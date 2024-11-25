#!/usr/bin/env python

from abc import ABC, abstractmethod
import math

#!/usr/bin/env python

class Region(ABC):
    @abstractmethod
    def contains(self, point):
        """
        Check if the given point (x, y, z) is inside the region.
        
        Parameters:
        - point: A tuple representing a point in space (x, y, z).
        
        Returns:
        - bool: True if the point is inside the region, False otherwise.
        """
        pass

def get_region(region_type: str, *args, **kwargs) -> Region:
    """
    Create an instance of a Region subclass based on the given type name.

    Parameters:
    - region_type: The name of the Region subclass to instantiate ('sphere' or 'block').
    - *args: Positional arguments to pass to the Region subclass constructor.
    - **kwargs: Keyword arguments to pass to the Region subclass constructor.

    Returns:
    - An instance of the specified Region subclass.

    Raises:
    - ValueError: If the specified region_type is not recognized.
    """
    region_classes = {
        'sphere': Sphere,
        'block': Block
    }

    if region_type not in region_classes:
        raise ValueError(f"Unknown region type: {region_type}")

    return region_classes[region_type](*args, **kwargs)

class Sphere(Region):
    def __init__(self, radius):
        """
        Initialize a Sphere instance with its center at (radius, radius, radius),
        making the x, y, and z planes tangent to the sphere.
        
        Parameters:
        - radius: The radius of the sphere. Must be a positive number.
        """
        assert radius > 0, "Radius must be a positive number"
        self.radius = radius
        # Set the center of the sphere to (radius, radius, radius)
        self.center = (radius, radius, radius)
        self.corner = (2*radius, 2*radius, 2*radius)

    def contains(self, point):
        """
        Check if the given point is inside the sphere.
        """
        dx = point[0] - self.center[0]
        dy = point[1] - self.center[1]
        dz = point[2] - self.center[2]
        distance_squared = dx**2 + dy**2 + dz**2
        return distance_squared < self.radius**2  # Points on the surface are not inside

class Block(Region):
    def __init__(self, corner):
        """
        Initialize a Block instance with one corner at the origin (0, 0, 0)
        and the opposite corner in the first quadrant.
        
        Parameters:
        - corner: A tuple representing the opposite corner of the block (x, y, z)
                  from the origin, ensuring x, y, z are positive.
        """
        assert all(c > 0 for c in corner), "Corner must be in the first quadrant"
        self.corner = corner

    def contains(self, point):
        """
        Check if the given point is inside the block defined by the origin
        and the corner in the first quadrant.
        """
        return all(0 < point[i] < self.corner[i] for i in range(3))

import os
import sys

TEST_MODE = (v:=os.environ.get('SIM_TEST_MODE','').upper()) \
    and v[0] in ('T', 'Y', '1')

if __name__ == '__main__':
    sys.exit(0) ### avoid running tests in production

if TEST_MODE:
    import pytest
    
    def test_block_contains_point_inside():
        block = Block((10, 10, 10))
        assert block.contains((5, 5, 5)) == True, \
        "The point is inside the block but was not recognized as such."
    
    def test_block_contains_point_outside():
        block = Block((10, 10, 10))
        assert block.contains((11, 11, 11)) == False, \
        "The point is outside the block but was incorrectly recognized as inside."
    
    def test_block_contains_point_on_edge():
        block = Block((10, 10, 10))
        assert block.contains((10, 10, 10)) == False, \
        "The point is on the edge of the block but was recognized as inside."
    
    def test_block_contains_point_at_origin():
        block = Block((10, 10, 10))
        assert block.contains((0, 0, 0)) == False, \
        "The origin point is outside the block but was not recognized as such."
    
    def test_block_corner_in_first_quadrant():
        with pytest.raises(AssertionError):
            Block((-1, -1, -1))  # This should raise an AssertionError because the corner is not in the first quadrant
    
    def test_block_contains_point_with_floats():
        block = Block((10.5, 10.5, 10.5))
        assert block.contains((10.1, 10.1, 10.1)) == True, \
        "The point is inside the block but was not recognized as such with floating point coordinates."

    def test_sphere_contains_point_inside():
        sphere = Sphere(5)
        # A point inside the sphere, not on the surface
        assert sphere.contains((6, 6, 4)) == True, \
        "The point is inside the sphere but was not recognized as such."
    
    def test_sphere_contains_point_outside():
        sphere = Sphere(5)
        # A point outside the sphere
        assert sphere.contains((11, 11, 11)) == False, \
        "The point is outside the sphere but was incorrectly recognized as inside."
    
    def test_sphere_contains_point_on_surface():
        sphere = Sphere(5)
        # A point exactly on the surface should not be considered inside
        assert sphere.contains((5, 5, 10)) == False, \
        "The point is on the surface of the sphere but was incorrectly recognized as inside."
    
    def test_sphere_center_correctly_set():
        sphere = Sphere(5)
        # The center should be at (radius, radius, radius)
        assert sphere.center == (5, 5, 5), \
        "The sphere's center is not correctly set at (radius, radius, radius)."
    
    def test_sphere_with_zero_radius():
        # Testing sphere initialization with zero radius should raise an AssertionError
        with pytest.raises(AssertionError):
            Sphere(0)
    
    def test_sphere_contains_point_near_origin():
        sphere = Sphere(5)
        # A point near the origin but outside the sphere
        assert sphere.contains((1, 1, 1)) == False, \
        "The point near the origin is incorrectly recognized as inside the sphere."

    def test_get_region_sphere():
        sphere = get_region('sphere', 5)
        assert isinstance(sphere, Sphere), "The created object is not an instance of Sphere"
        assert sphere.radius == 5, "The radius of the created sphere is incorrect"

    def test_get_region_block():
        block = get_region('block', (10, 10, 10))
        assert isinstance(block, Block), "The created object is not an instance of Block"
        assert block.corner == (10, 10, 10), "The corner of the created block is incorrect"

    def test_get_region_invalid_type():
        with pytest.raises(ValueError, match="Unknown region type: invalid"):
            get_region('invalid', 5)

    def test_get_region_block_with_floats():
        block = get_region('block', (10.5, 10.5, 10.5))
        assert isinstance(block, Block), "The created object is not an instance of Block"
        assert block.corner == (10.5, 10.5, 10.5), "The corner of the created block is incorrect"

    def test_get_region_sphere_with_zero_radius():
        with pytest.raises(AssertionError):
            get_region('sphere', 0)

    def test_get_region_block_with_negative_corner():
        with pytest.raises(AssertionError):
            get_region('block', (-1, -1, -1))