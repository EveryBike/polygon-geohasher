import unittest
import time
from shapely import to_geojson
from polygon_geohasher.polygon_geohasher import (
    polygon_to_geohashes,
    geohash_to_polygon,
    geohashes_to_polygon,
    polygon_to_geohashes_with_intersection,
)
from shapely import geometry
from shapely.ops import unary_union


class TestBasicMethods(unittest.TestCase):
    def test_one_geohash(self):
        test_geohash = "x1"
        test_polygon = geohash_to_polygon(test_geohash)
        polygon = geohashes_to_polygon(polygon_to_geohashes(test_polygon, 2))

        self.assertEqual(test_polygon.area, polygon.area)
        self.assertTrue(test_polygon.covers(test_polygon.intersection(polygon)))
        self.assertTrue(test_polygon.intersection(polygon).covers(test_polygon))

    def test_triangle(self):
        test_polygon = geometry.Polygon(
            [
                (-99.1795917, 19.432134),
                (-99.1656847, 19.429034),
                (-99.1776492, 19.414236),
                (-99.1795917, 19.432134),
            ]
        )

        polygon = geohashes_to_polygon(polygon_to_geohashes(test_polygon, 7))
        self.assertTrue(polygon.area <= test_polygon.area)
        self.assertTrue(
            all(polygon.covers(geometry.Point(c)) for c in polygon.boundary.coords)
        )

        polygon = geohashes_to_polygon(polygon_to_geohashes(test_polygon, 7, False))
        self.assertTrue(polygon.area >= test_polygon.area)


class TestIntersectionCorrectness(unittest.TestCase):
    """Verify that the union of intersected cells faithfully reproduces the input polygon."""

    def _assert_union_matches_polygon(self, polygon, precision, grid_size=0, tolerance=1e-8):
        """Helper: union of all intersection values should match the original polygon area."""
        result = polygon_to_geohashes_with_intersection(polygon, precision, False, grid_size)
        self.assertTrue(len(result) > 0, "Expected at least one geohash")

        parts = [v for v in result.values() if not v.is_empty]
        union = unary_union(parts)
        self.assertAlmostEqual(
            polygon.area, union.area, places=7,
            msg=f"Area mismatch: polygon={polygon.area}, union={union.area}"
        )
        # The union should cover the original polygon within tolerance
        diff = polygon.difference(union)
        self.assertLess(diff.area, tolerance, f"Uncovered area: {diff.area}")

    def test_rectangle(self):
        """Axis-aligned rectangle — simplest case."""
        poly = geometry.box(127.01, 37.39, 127.02, 37.40)
        self._assert_union_matches_polygon(poly, 7)

    def test_triangle(self):
        """Non-rectangular polygon with diagonal edges."""
        poly = geometry.Polygon([
            (-99.1795917, 19.432134),
            (-99.1656847, 19.429034),
            (-99.1776492, 19.414236),
            (-99.1795917, 19.432134),
        ])
        self._assert_union_matches_polygon(poly, 7, grid_size=0.000001, tolerance=0.000001)

    def test_l_shaped_polygon(self):
        """Concave L-shape — tests correct handling of non-convex geometry."""
        poly = geometry.Polygon([
            (127.00, 37.40),
            (127.02, 37.40),
            (127.02, 37.42),
            (127.01, 37.42),
            (127.01, 37.41),
            (127.00, 37.41),
            (127.00, 37.40),
        ])
        self._assert_union_matches_polygon(poly, 6)

    def test_narrow_strip(self):
        """Very narrow polygon — many boundary cells, few fully contained."""
        poly = geometry.box(127.00, 37.400, 127.10, 37.401)
        self._assert_union_matches_polygon(poly, 7)

    def test_small_polygon_high_precision(self):
        """Polygon smaller than a single geohash cell at low precision."""
        poly = geometry.box(127.010, 37.390, 127.011, 37.391)
        self._assert_union_matches_polygon(poly, 8)

    def test_precision_levels(self):
        """Same polygon at different precisions should all produce correct area."""
        poly = geometry.box(127.01, 37.39, 127.03, 37.41)
        for precision in [4, 5, 6, 7]:
            with self.subTest(precision=precision):
                self._assert_union_matches_polygon(poly, precision)


class TestIntersectionInnerMode(unittest.TestCase):
    """Verify inner=True: returned cells must be fully contained within the polygon."""

    def test_inner_cells_are_contained(self):
        poly = geometry.Polygon([
            (-99.1795917, 19.432134),
            (-99.1656847, 19.429034),
            (-99.1776492, 19.414236),
            (-99.1795917, 19.432134),
        ])
        result = polygon_to_geohashes_with_intersection(poly, 7, inner=True)
        for gh, cell_poly in result.items():
            self.assertTrue(
                poly.contains(cell_poly),
                f"Geohash {gh} is not fully contained in polygon"
            )

    def test_inner_subset_of_outer(self):
        """inner=True keys must be a subset of inner=False keys."""
        poly = geometry.box(127.01, 37.39, 127.02, 37.40)
        inner = polygon_to_geohashes_with_intersection(poly, 7, inner=True)
        outer = polygon_to_geohashes_with_intersection(poly, 7, inner=False)
        self.assertTrue(
            set(inner.keys()).issubset(set(outer.keys())),
            "Inner geohashes should be a subset of outer geohashes"
        )

    def test_inner_area_less_than_polygon(self):
        """Union of inner cells must not exceed the polygon area."""
        poly = geometry.Polygon([
            (-99.18, 19.43),
            (-99.16, 19.43),
            (-99.17, 19.41),
            (-99.18, 19.43),
        ])
        result = polygon_to_geohashes_with_intersection(poly, 7, inner=True)
        if len(result) > 0:
            union = unary_union(list(result.values()))
            self.assertLessEqual(union.area, poly.area + 1e-12)


class TestIntersectionBoundaryCells(unittest.TestCase):
    """Verify that boundary cells have geometrically correct intersections."""

    def test_boundary_cells_within_polygon(self):
        """Each intersection polygon should be contained within the original polygon."""
        poly = geometry.box(127.01, 37.39, 127.02, 37.40)
        result = polygon_to_geohashes_with_intersection(poly, 7, inner=False, grid_size=0)
        for gh, intersection in result.items():
            self.assertTrue(
                poly.buffer(1e-9).contains(intersection),
                f"Intersection of {gh} extends outside the polygon"
            )

    def test_boundary_cells_within_geohash(self):
        """Each intersection polygon should be contained within its geohash cell."""
        poly = geometry.box(127.01, 37.39, 127.02, 37.40)
        result = polygon_to_geohashes_with_intersection(poly, 7, inner=False, grid_size=0)
        for gh, intersection in result.items():
            cell = geohash_to_polygon(gh)
            self.assertTrue(
                cell.buffer(1e-9).contains(intersection),
                f"Intersection of {gh} extends outside its geohash cell"
            )

    def test_fully_contained_cells_equal_geohash(self):
        """Cells fully inside the polygon should return the full geohash polygon."""
        poly = geometry.box(127.00, 37.38, 127.04, 37.42)
        result = polygon_to_geohashes_with_intersection(poly, 7, inner=False, grid_size=0)
        for gh, intersection in result.items():
            cell = geohash_to_polygon(gh)
            if poly.contains(cell):
                self.assertAlmostEqual(
                    intersection.area, cell.area, places=12,
                    msg=f"Fully contained cell {gh} should return full cell area"
                )


class TestIntersectionEdgeCases(unittest.TestCase):

    def test_empty_polygon(self):
        poly = geometry.Polygon()
        result = polygon_to_geohashes_with_intersection(poly, 7, inner=False)
        self.assertEqual(len(result), 0)

    def test_grid_size_none(self):
        """grid_size=None should work without error."""
        poly = geometry.box(127.01, 37.39, 127.02, 37.40)
        result = polygon_to_geohashes_with_intersection(poly, 7, inner=False, grid_size=None)
        self.assertTrue(len(result) > 0)

    def test_grid_size_nonzero(self):
        """Non-zero grid_size should snap coordinates and still produce valid results."""
        poly = geometry.Polygon([
            (-99.1795917, 19.432134),
            (-99.1656847, 19.429034),
            (-99.1776492, 19.414236),
            (-99.1795917, 19.432134),
        ])
        result = polygon_to_geohashes_with_intersection(poly, 7, inner=False, grid_size=0.000001)
        parts = [v for v in result.values() if not v.is_empty]
        union = unary_union(parts)
        self.assertAlmostEqual(poly.area, union.area, places=5)

    def test_no_overlapping_intersections(self):
        """Intersection polygons should not significantly overlap each other."""
        poly = geometry.box(127.01, 37.39, 127.02, 37.40)
        result = polygon_to_geohashes_with_intersection(poly, 6, inner=False, grid_size=0)
        values = [v for v in result.values() if not v.is_empty]
        sum_of_areas = sum(v.area for v in values)
        union_area = unary_union(values).area
        self.assertAlmostEqual(
            sum_of_areas, union_area, places=10,
            msg=f"Overlapping intersections detected: sum={sum_of_areas}, union={union_area}"
        )


class TestBenchmark(unittest.TestCase):
    def test_benchmark(self):
        """Hierarchical approach should handle large polygons efficiently."""
        poly = geometry.Polygon([
            (-99.2, 19.3),
            (-99.1, 19.3),
            (-99.1, 19.5),
            (-99.2, 19.5),
            (-99.2, 19.3),
        ])

        start = time.time()
        result = polygon_to_geohashes_with_intersection(poly, 7, False, 0)
        elapsed = time.time() - start

        print(f"\nBenchmark: {len(result)} geohashes in {elapsed:.3f}s")
        self.assertTrue(len(result) > 0)
        self.assertTrue(elapsed < 30, f"Too slow: {elapsed:.3f}s")


if __name__ == "__main__":
    unittest.main()
