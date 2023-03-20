import unittest
from functools import reduce
from shapely import to_geojson
from polygon_geohasher.polygon_geohasher import (
    polygon_to_geohashes,
    geohash_to_polygon,
    geohashes_to_polygon,
    polygon_to_geohashes_with_intersection,
)
from shapely import geometry


class TestSimpleMethods(unittest.TestCase):
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
        
        grid_size = 0.000001
        intersections = polygon_to_geohashes_with_intersection(test_polygon, 7, False, grid_size)
        multi_polygon = geometry.MultiPolygon(intersections.values())
        print(to_geojson(test_polygon))
        print(to_geojson(multi_polygon))
        # check visually on https://geojson.io/
        self.assertTrue(test_polygon.area - multi_polygon.area <= grid_size)
        # print(test_polygon.area)
        # print(multi_polygon.area)


if __name__ == "__main__":
    unittest.main()
