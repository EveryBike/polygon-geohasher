import geohash
import queue

from shapely import geometry
from shapely.ops import unary_union


def geohash_to_polygon(geo):
    """
    :param geo: String that represents the geohash.
    :return: Returns a Shapely's Polygon instance that represents the geohash.
    """
    lat_centroid, lng_centroid, lat_offset, lng_offset = geohash.decode_exactly(geo)

    corner_1 = (lat_centroid - lat_offset, lng_centroid - lng_offset)[::-1]
    corner_2 = (lat_centroid - lat_offset, lng_centroid + lng_offset)[::-1]
    corner_3 = (lat_centroid + lat_offset, lng_centroid + lng_offset)[::-1]
    corner_4 = (lat_centroid + lat_offset, lng_centroid - lng_offset)[::-1]

    return geometry.Polygon([corner_1, corner_2, corner_3, corner_4, corner_1])


def polygon_to_geohashes(polygon, precision, inner=True):
    """
    :param polygon: shapely polygon.
    :param precision: int. Geohashes' precision that form resulting polygon.
    :param inner: bool, default 'True'. If false, geohashes that are completely outside from the polygon are ignored.
    :return: set. Set of geohashes that form the polygon.
    """
    inner_geohashes = set()
    outer_geohashes = set()

    envelope = polygon.envelope
    centroid = polygon.centroid

    testing_geohashes = queue.Queue()
    testing_geohashes.put(geohash.encode(centroid.y, centroid.x, precision))

    while not testing_geohashes.empty():
        current_geohash = testing_geohashes.get()

        if (
            current_geohash not in inner_geohashes
            and current_geohash not in outer_geohashes
        ):
            current_polygon = geohash_to_polygon(current_geohash)

            condition = (
                envelope.contains(current_polygon)
                if inner
                else envelope.intersects(current_polygon)
            )

            if condition:
                if inner:
                    if polygon.contains(current_polygon):
                        inner_geohashes.add(current_geohash)
                    else:
                        outer_geohashes.add(current_geohash)
                else:
                    if polygon.intersects(current_polygon):
                        inner_geohashes.add(current_geohash)
                    else:
                        outer_geohashes.add(current_geohash)
                for neighbor in geohash.neighbors(current_geohash):
                    if (
                        neighbor not in inner_geohashes
                        and neighbor not in outer_geohashes
                    ):
                        testing_geohashes.put(neighbor)

    return inner_geohashes

def polygon_to_geohashes_with_intersection(polygon, precision, inner=True, grid_size=None):
    """
    :param polygon: shapely polygon.
    :param precision: int. Geohashes' precision that form resulting polygon.
    :param inner: bool, default 'True'. If false, geohashes that are completely outside from the polygon are ignored.
    :param grid_size: double, default 0. Precision grid size, 
    if none-zero, nput coordinates will be snapped to a precision grid of that size and resulting coordinates will be snapped to that same grid
    If 0, this operation will use double precision coordinates. 
    If None, the highest precision of the inputs will be used.
    :return: dict. key is the geohash value is the intersected polygon.
    """
    inner_geohashes = {}
    outer_geohashes = set()
    
    envelope = polygon.envelope
    centroid = polygon.centroid

    testing_geohashes = queue.Queue()
    testing_geohashes.put(geohash.encode(centroid.y, centroid.x, precision))

    while not testing_geohashes.empty():
        current_geohash = testing_geohashes.get()

        if (
            current_geohash not in inner_geohashes
            and current_geohash not in outer_geohashes
        ):
            current_polygon = geohash_to_polygon(current_geohash)

            condition = (
                envelope.contains(current_polygon)
                if inner
                else envelope.intersects(current_polygon)
            )

            if condition:
                if inner:
                    if polygon.contains(current_polygon):
                        inner_geohashes[current_geohash] = current_polygon
                    else:
                        outer_geohashes.add(current_geohash)
                else:
                    if polygon.intersects(current_polygon):
                        inner_geohashes[current_geohash] = polygon.intersection(current_polygon, grid_size=grid_size)
                    else:
                        outer_geohashes.add(current_geohash)
                for neighbor in geohash.neighbors(current_geohash):
                    if (
                        neighbor not in inner_geohashes
                        and neighbor not in outer_geohashes
                    ):
                        testing_geohashes.put(neighbor)

    return inner_geohashes

def geohashes_to_polygon(geohashes):
    """
    :param geohashes: array-like. List of geohashes to form resulting polygon.
    :return: shapely geometry. Resulting Polygon after combining geohashes.
    """
    return unary_union([geohash_to_polygon(g) for g in geohashes])
