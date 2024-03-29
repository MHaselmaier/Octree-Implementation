package de.fhkl.imst.i.cgma.raytracer.octree;

import java.awt.Color;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.Vector;

import de.fhkl.imst.i.cgma.raytracer.file.I_Sphere;
import de.fhkl.imst.i.cgma.raytracer.file.RT_Object;
import de.fhkl.imst.i.cgma.raytracer.file.T_Mesh;

public class Octree
{
	private static final double EPSILON = 1E-5;
	private static final int MINIMUM_NUMBER_OF_TRIANGLES = 30;
	private static final float[] VOXEL_RIGHT_FACE_NORMAL = {1, 0, 0};
	private static final float[] VOXEL_LEFT_FACE_NORMAL = {-1, 0, 0};
	private static final float[] VOXEL_TOP_FACE_NORMAL = {0, 1, 0};
	private static final float[] VOXEL_BOTTOM_FACE_NORMAL = {0, -1, 0};
	private static final float[] VOXEL_FRONT_FACE_NORMAL = {0, 0, 1};
	private static final float[] VOXEL_BACK_FACE_NORMAL = {0, 0, -1};

	private float[] ambientLightColor = { 0.25f, 0.25f, 0.25f };
	private float[] specularLight = { 1.0f, 1.0f, 1.0f };
	private float[] centerOfPointLight = { 4.0f, 4.0f, 2.0f };
	
	private float[] min = new float[3];
	private float[] max = new float[3];
	private Octree parent;
	private Map<RT_Object, Integer[]> objects;
	private Octree[] children;

	public Octree(Vector<RT_Object> objects)
	{
		Octree.calculateBoundingBox(objects, this.min, this.max);
		this.objects = Octree.generateObjectTriangleIndexMap(objects);
		this.children = generateChildren();
	}
	
	private static void calculateBoundingBox(Vector<RT_Object> objects, float[] min, float[] max)
    {
        Arrays.fill(min, Float.POSITIVE_INFINITY);
        Arrays.fill(max, Float.NEGATIVE_INFINITY);
        for (RT_Object object: objects)
        {
            for(int i = 0; 3 > i; ++i)
            {
                if(min[i] > object.min[i]) min[i] = object.min[i];
                if(max[i] < object.max[i]) max[i] = object.max[i];
            }
        }
    }
	
	private static Map<RT_Object, Integer[]> generateObjectTriangleIndexMap(Vector<RT_Object> objects)
    {
        Map<RT_Object, Integer[]> containedObjects = new HashMap<>();
        for (RT_Object object: objects)
        {
            Integer[] indices = null;
            if (object instanceof T_Mesh)
            {
                indices = new Integer[((T_Mesh)object).triangles.length];
            }
            else if(object instanceof I_Sphere)
            {
                indices = new Integer[1];
            }
            for (int i = 0; indices.length > i; ++i)
            {
                indices[i] = i;
            }
            containedObjects.put(object, indices);
        }
        return containedObjects;
    }
	
	private Octree[] generateChildren()
	{
		if (Octree.MINIMUM_NUMBER_OF_TRIANGLES >= this.countContainedTriangles()) return null;

		float[] delta = new float[3];
		delta[0] = (this.max[0] - this.min[0]) / 2f;
		delta[1] = (this.max[1] - this.min[1]) / 2f;
		delta[2] = (this.max[2] - this.min[2]) / 2f;

		Octree[] children = new Octree[8];
		float[] newMin;
		float[] newMax;

		// Bottom-Back-Left
		newMin = new float[] { this.min[0], this.min[1], this.min[2] };
		newMax = new float[] { this.min[0] + delta[0], this.min[1] + delta[1], this.min[2] + delta[2] };
		children[0] = new Octree(this, newMin, newMax);

		// Bottom-Back-Right
		newMin = new float[] { this.min[0] + delta[0], this.min[1], this.min[2] };
		newMax = new float[] { this.max[0], this.min[1] + delta[1], this.min[2] + delta[2] };
		children[1] = new Octree(this, newMin, newMax);

		// Bottom-Front-Right
		newMin = new float[] { this.min[0] + delta[0], this.min[1], this.min[2] + delta[2] };
		newMax = new float[] { this.max[0], this.min[1] + delta[1], this.max[2] };
		children[2] = new Octree(this, newMin, newMax);

		// Bottom-Front-Left
		newMin = new float[] { this.min[0], this.min[1], this.min[2] + delta[2] };
		newMax = new float[] { this.min[0] + delta[0], this.min[1] + delta[1], this.max[2] };
		children[3] = new Octree(this, newMin, newMax);

		// Top-Back-Left
		newMin = new float[] { this.min[0], this.min[1] + delta[1], this.min[2] };
		newMax = new float[] { this.min[0] + delta[0], this.max[1], this.min[2] + delta[2] };
		children[4] = new Octree(this, newMin, newMax);

		// Top-Back-Right
		newMin = new float[] { this.min[0] + delta[0], this.min[1] + delta[1], this.min[2] };
		newMax = new float[] { this.max[0], this.max[1], this.min[2] + delta[2] };
		children[5] = new Octree(this, newMin, newMax);

		// Top-Front-Right
		newMin = new float[] { this.min[0] + delta[0], this.min[1] + delta[1], this.min[2] + delta[2] };
		newMax = new float[] { this.max[0], this.max[1], this.max[2] };
		children[6] = new Octree(this, newMin, newMax);

		// Top-Front-Left
		newMin = new float[] { this.min[0], this.min[1] + delta[1], this.min[2] + delta[2] };
		newMax = new float[] { this.min[0] + delta[0], this.max[1], this.max[2] };
		children[7] = new Octree(this, newMin, newMax);
		
		return children;
	}
	
	private int countContainedTriangles()
	{
		int containedTriangles = 0;
		for (Integer[] indices : this.objects.values())
		{
			containedTriangles += indices.length;
		}
		return containedTriangles;
	}

	private Octree(Octree parent, float[] min, float[] max)
	{
		this.min = min;
		this.max = max;
		this.parent = parent;
		this.objects = findContainedObjects();
		this.children = generateChildren();
	}

	private Map<RT_Object, Integer[]> findContainedObjects()
	{
		Map<RT_Object, Integer[]> containedObjects = new HashMap<>();

		for (Map.Entry<RT_Object, Integer[]> object : this.parent.objects.entrySet())
		{
			RT_Object scene = object.getKey();

			if (!couldObjectBeInsideVoxel(scene)) continue;

			if (scene instanceof I_Sphere)
			{
				// Nur Bounding Box Check ...
				// Ist ja eh implizit und besteht nicht aus Dreiecken.
				containedObjects.put(scene, new Integer[1]);
			}
			else if (scene instanceof T_Mesh)
			{
				T_Mesh mesh = (T_Mesh) scene;
				Set<Integer> indices = new TreeSet<>();
				// Teste für jedes Dreieck des Objektes,
				// ob es innerhalb des Voxels liegt.
				for (int i = 0; mesh.triangles.length > i; ++i)
				{
					if (isTriangleContained(new Triangle(mesh, i)))
					{
						indices.add(i);
					}
				}
				if (0 < indices.size())
				{
					containedObjects.put(mesh, indices.toArray(new Integer[0]));
				}
			}
		}
		return containedObjects;
	}
	
	private boolean couldObjectBeInsideVoxel(RT_Object scene) {
		return (Octree.areBoxPointsInsideBox(scene.min, scene.max, this.min, this.max) ||
				Octree.areBoxPointsInsideBox(this.min, this.max, scene.min, scene.max) ||
				Octree.areBoxesIntersecting(scene.min, scene.max, this.min, this.max) ||
				Octree.areBoxesIntersecting(this.min, this.max, scene.min, scene.max));
	}
	
	private static boolean areBoxPointsInsideBox(float[] firstMin, float[] firstMax, float[] secondMin, float[] secondMax)
	{
		return (Octree.isPointInsideBox(firstMin[0], firstMin[1], firstMin[2], secondMin, secondMax) ||
				Octree.isPointInsideBox(firstMin[0], firstMin[1], firstMax[2], secondMin, secondMax) ||
				Octree.isPointInsideBox(firstMin[0], firstMax[1], firstMin[2], secondMin, secondMax) ||
				Octree.isPointInsideBox(firstMin[0], firstMax[1], firstMax[2], secondMin, secondMax) ||
				Octree.isPointInsideBox(firstMax[0], firstMin[1], firstMin[2], secondMin, secondMax) ||
				Octree.isPointInsideBox(firstMax[0], firstMin[1], firstMax[2], secondMin, secondMax) ||
				Octree.isPointInsideBox(firstMax[0], firstMax[1], firstMin[2], secondMin, secondMax) ||
				Octree.isPointInsideBox(firstMax[0], firstMax[1], firstMax[2], secondMin, secondMax));
	}
	
	private static boolean isPointInsideBox(float x, float y, float z, float[] min, float[] max)
	{
		return (((min[0] < x || EPSILON > Math.abs(min[0] - x)) && (max[0] > x || EPSILON > Math.abs(max[0] - x))) &&
				((min[1] < y || EPSILON > Math.abs(min[1] - y)) && (max[1] > y || EPSILON > Math.abs(max[1] - y))) &&
				((min[2] < z || EPSILON > Math.abs(min[2] - z)) && (max[2] > z || EPSILON > Math.abs(max[2] - z))));
	}
	
	private static boolean areBoxesIntersecting(float[] firstMin, float[] firstMax, float[] secondMin, float[] secondMax)
	{
		return (((firstMin[0] <= secondMin[0] && firstMax[0] >= secondMax[0]) &&
				 (firstMin[1] >= secondMin[1] && firstMax[1] <= secondMax[1]) &&
				 (firstMin[2] >= secondMin[2] && firstMax[2] <= secondMax[2])) ||
				((firstMin[0] >= secondMin[0] && firstMax[0] <= secondMax[0]) &&
				 (firstMin[1] <= secondMin[1] && firstMax[1] >= secondMax[1]) &&
				 (firstMin[2] >= secondMin[2] && firstMax[2] <= secondMax[2])) ||
				((firstMin[0] >= secondMin[0] && firstMax[0] <= secondMax[0]) &&
				 (firstMin[1] >= secondMin[1] && firstMax[1] <= secondMax[1]) &&
				 (firstMin[2] <= secondMin[2] && firstMax[2] >= secondMax[2])));
	}
	
	private boolean isTriangleContained(Triangle triangle)
	{
		return (Octree.isPointInsideBox(triangle.firstPoint, this.min, this.max) ||
				Octree.isPointInsideBox(triangle.secondPoint, this.min, this.max) ||
				Octree.isPointInsideBox(triangle.thirdPoint, this.min, this.max) ||
				isTriangleEdgeIntersectingVoxel(triangle) ||
				isVoxelEdgeIntersectingTriangle(triangle));
	}
	
	private static boolean isPointInsideBox(float[] point, float[] min, float[] max)
	{
		return Octree.isPointInsideBox(point[0], point[1], point[2], min, max);
	}
	
	private boolean isTriangleEdgeIntersectingVoxel(Triangle triangle)
	{
		float[] p1p2 = { triangle.secondPoint[0] - triangle.firstPoint[0],
						 triangle.secondPoint[1] - triangle.firstPoint[1],
						 triangle.secondPoint[2] - triangle.firstPoint[2] };
		float[] p1p3 = { triangle.thirdPoint[0] - triangle.firstPoint[0],
						 triangle.thirdPoint[1] - triangle.firstPoint[1],
						 triangle.thirdPoint[2] - triangle.firstPoint[2] };
		float[] p2p3 = { triangle.thirdPoint[0] - triangle.secondPoint[0],
						 triangle.thirdPoint[1] - triangle.secondPoint[1],
						 triangle.thirdPoint[2] - triangle.secondPoint[2] };

		return (isTriangleEdgeIntersectingFace(triangle.firstPoint, triangle.secondPoint, p1p2, p1p3, p2p3, Octree.VOXEL_RIGHT_FACE_NORMAL, this.max) ||
				isTriangleEdgeIntersectingFace(triangle.firstPoint, triangle.secondPoint, p1p2, p1p3, p2p3, Octree.VOXEL_LEFT_FACE_NORMAL, this.min) ||
				isTriangleEdgeIntersectingFace(triangle.firstPoint, triangle.secondPoint, p1p2, p1p3, p2p3, Octree.VOXEL_TOP_FACE_NORMAL, this.max) ||
				isTriangleEdgeIntersectingFace(triangle.firstPoint, triangle.secondPoint, p1p2, p1p3, p2p3, Octree.VOXEL_BOTTOM_FACE_NORMAL, this.min) ||
				isTriangleEdgeIntersectingFace(triangle.firstPoint, triangle.secondPoint, p1p2, p1p3, p2p3, Octree.VOXEL_FRONT_FACE_NORMAL, this.max) ||
				isTriangleEdgeIntersectingFace(triangle.firstPoint, triangle.secondPoint, p1p2, p1p3, p2p3, Octree.VOXEL_BACK_FACE_NORMAL, this.min));
	}
	
	private boolean isTriangleEdgeIntersectingFace(float[] p1, float[] p2, float[] p1p2, float[] p1p3, float[] p2p3, float[] n, float[] facePoint)
	{
		float[] intersectionPoint = new float[3];
		return ((isIntersectingInRange(p1, p1p2, n, facePoint, intersectionPoint) && Octree.isPointInsideBox(intersectionPoint, this.min, this.max)) ||
				(isIntersectingInRange(p1, p1p3, n, facePoint, intersectionPoint) && Octree.isPointInsideBox(intersectionPoint, this.min, this.max)) ||
				(isIntersectingInRange(p2, p2p3, n, facePoint, intersectionPoint) && Octree.isPointInsideBox(intersectionPoint, this.min, this.max)));
	}
	
	private static boolean isIntersectingInRange(float[] start, float[] direction, float[] normal, float[] point, float[] intersectionPoint)
	{
		if (Octree.isIntersecting(start, direction, normal, point, intersectionPoint))
		{
			float factor = Octree.calculateIntersectionPointFactor(start, direction, intersectionPoint);
			return (1 >= factor && 0 <= factor);
		}

		return false;
	}
	
	private static float calculateIntersectionPointFactor(float[] start, float[] direction, float[] intersectionPoint)
	{
		if (0 < Math.abs(direction[0])) return (intersectionPoint[0] - start[0]) / direction[0];
		if (0 < Math.abs(direction[1])) return (intersectionPoint[1] - start[1]) / direction[1];
		if (0 < Math.abs(direction[2])) return (intersectionPoint[2] - start[2]) / direction[2];
		return Float.NaN;
	}
	
	private static boolean isIntersecting(float[] start, float[] direction, float[] normal, float[] point, float[] intersectionPoint)
	{
		float vn = direction[0] * normal[0] + direction[1] * normal[1] + direction[2] * normal[2];
		if (Math.abs(vn) < EPSILON) return false;

		float pen = (point[0] - start[0]) * normal[0] + (point[1] - start[1]) * normal[1] + (point[2] - start[2]) * normal[2];
		float t = pen / vn;

		intersectionPoint[0] = start[0] + t * direction[0];
		intersectionPoint[1] = start[1] + t * direction[1];
		intersectionPoint[2] = start[2] + t * direction[2];

		return true;
	}
	
	private boolean isVoxelEdgeIntersectingTriangle(Triangle triangle)
	{
		float[] bottomBackLeft = { this.min[0], this.min[1], this.min[2] };
		float[] bottomFrontRight = { this.max[0], this.min[1], this.max[2] };
		float[] topBackRight = { this.max[0], this.max[1], this.min[2] };
		float[] topFrontLeft = { this.min[0], this.max[1], this.max[2] };

		return (isIntersectingTriangle(bottomBackLeft, new float[] { this.max[0] - this.min[0], 0, 0 }, triangle) ||
				isIntersectingTriangle(bottomBackLeft, new float[] { 0, this.max[1] - this.min[1], 0 }, triangle) ||
				isIntersectingTriangle(bottomBackLeft, new float[] { 0, 0, this.max[2] - this.min[2] }, triangle) ||
				isIntersectingTriangle(bottomFrontRight, new float[] { this.min[0] - this.max[0], 0, 0 }, triangle) ||
				isIntersectingTriangle(bottomFrontRight, new float[] { 0, this.max[1] - this.min[1], 0 }, triangle) ||
				isIntersectingTriangle(bottomFrontRight, new float[] { 0, 0, this.min[2] - this.max[2] }, triangle) ||
				isIntersectingTriangle(topBackRight, new float[] { this.min[0] - this.max[0], 0, 0 }, triangle) ||
				isIntersectingTriangle(topBackRight, new float[] { 0, this.min[1] - this.max[1], 0 }, triangle) ||
				isIntersectingTriangle(topBackRight, new float[] { 0, 0, this.max[2] - this.min[2] }, triangle) ||
				isIntersectingTriangle(topFrontLeft, new float[] { this.max[0] - this.min[0], 0, 0 }, triangle) ||
				isIntersectingTriangle(topFrontLeft, new float[] { 0, this.min[1] - this.max[1], 0 }, triangle) ||
				isIntersectingTriangle(topFrontLeft, new float[] { 0, 0, this.min[2] - this.max[2] }, triangle));
	}
	
	private boolean isIntersectingTriangle(float[] start, float[] direction, Triangle triangle)
	{
		float[] intersectionPoint = new float[3];
		return (isIntersectingInRange(start, direction, triangle.normal, triangle.firstPoint, intersectionPoint) &&
				triangleTest(intersectionPoint, triangle));
	}
	
	private boolean triangleTest(float[] point, Triangle triangle)
	{
		return triangleTest(point, triangle.firstPoint, triangle.secondPoint, triangle.thirdPoint, triangle.area, new float[3]);
	}
	
	private boolean triangleTest(float[] point, float[] p1, float[] p2, float[] p3, float a, float ai[])
	{
		ai[0] = calculateArea(p1, p2, point);
		ai[1] = calculateArea(point, p2, p3);
		ai[2] = calculateArea(p1, point, p3);

		return (EPSILON > Math.abs(a - (ai[0] + ai[1] + ai[2])));
	}
	
	private float calculateArea(float[] p1, float[] p2, float[] p3)
	{
		// a = Vi2-Vi1, b = Vi3-Vi1
		float ax, ay, az, bx, by, bz;
		ax = p2[0] - p1[0];
		ay = p2[1] - p1[1];
		az = p2[2] - p1[2];

		bx = p3[0] - p1[0];
		by = p3[1] - p1[1];
		bz = p3[2] - p1[2];

		// n = a x b
		float[] n = new float[3];
		n[0] = ay * bz - az * by;
		n[1] = -(ax * bz - az * bx);
		n[2] = ax * by - ay * bx;

		// normalize n, calculate and return area of triangle
		return normalize(n) / 2;
	}
	
	private float normalize(float[] vector)
	{
		float length = (float)Math.sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
		vector[0] /= length;
		vector[1] /= length;
		vector[2] /= length;
		return length;
	}
	
	public Color traceRay(float[] eye, float[] ray)
	{
		if (null == this.children) return traceRayInsideVoxel(eye, ray);
		
		Color color = null;
		float smallestFactor = 0;
		do
		{
			Octree voxel = null;
			float nextFactor = Float.POSITIVE_INFINITY;
			for (int i = 0; this.children.length > i; ++i)
			{
				float factor = this.children[i].calculateIntersectionPointRayFactor(eye, ray);
				if (nextFactor > factor && smallestFactor < factor)
				{
					nextFactor = factor;
					voxel = this.children[i];
				}
			}
			smallestFactor = nextFactor;
			if (null == voxel)
			{
				return null;
			}
			color = voxel.traceRay(eye, ray);
		} while (null == color);
		return color;
	}

	private Color traceRayInsideVoxel(float[] eye, float[] ray)
	{
		if (null != this.children) return null;

		float minT = Float.MAX_VALUE;
		RT_Object minObjectsKey = null;
		int minIndex = -1;
		float[] ip = new float[3];
		float[] minIP = new float[3];
		float[] minN = new float[3];
		float[] minMaterial = new float[3];
		float minMaterialN = 1;
		float bu = 0, bv = 0, bw = 1;

		float[] v = new float[3];
		float[] l = new float[3];

		// viewing vector at intersection point
		v[0] = eye[0] - ray[0];
		v[1] = eye[1] - ray[1];
		v[2] = eye[2] - ray[2];
		normalize(v);

		RT_Object scene;
		I_Sphere sphere;
		T_Mesh mesh;

		for (Map.Entry<RT_Object, Integer[]> object : this.objects.entrySet())
		{
			scene = object.getKey();
			if (scene instanceof I_Sphere)
			{
				sphere = (I_Sphere) scene;

				// no bounding box hit? -> next object
				if (!Octree.boundingBoxHit(sphere, eye, ray)) continue;

				// ray intersection uses quadratic equation
				float a, b, c, d;
				a = ray[0] * ray[0] + ray[1] * ray[1] + ray[2] * ray[2];
				b = 2f * (ray[0] * (eye[0] - sphere.center[0]) + ray[1] * (eye[1] - sphere.center[1])
						+ ray[2] * (eye[2] - sphere.center[2]));
				c = ((eye[0] - sphere.center[0]) * (eye[0] - sphere.center[0])
						+ (eye[1] - sphere.center[1]) * (eye[1] - sphere.center[1])
						+ (eye[2] - sphere.center[2]) * (eye[2] - sphere.center[2]))
						- sphere.radius * sphere.radius;

				// positive discriminant determines intersection
				d = b * b - 4f * a * c;
				// no intersection point? => next object
				if (d <= 0) continue;

				// from here: intersection takes place!
				// calculate first intersection point with sphere along the ray
				float t = (float) (-b - Math.sqrt(d)) / (2 * a);

				// already a closer intersection point? => next object
				if (t >= minT) continue;

				ip[0] = eye[0] + t * ray[0];
				ip[1] = eye[1] + t * ray[1];
				ip[2] = eye[2] + t * ray[2];
				
				// ### CHANGED ###
				if (!Octree.isPointInsideBox(ip, this.min, this.max)) continue;

				// from here: t < minT
				// I'm the winner until now!
				minT = t;
				minObjectsKey = scene;

				// prepare everything for phong shading
				// the intersection point
				minIP[0] = ip[0];
				minIP[1] = ip[1];
				minIP[2] = ip[2];

				// the normal vector at the intersection point
				minN[0] = minIP[0] - sphere.center[0];
				minN[1] = minIP[1] - sphere.center[1];
				minN[2] = minIP[2] - sphere.center[2];
				normalize(minN);

				// the material
				minMaterial = sphere.material;
				minMaterialN = sphere.materialN;
			}
			else if (scene instanceof T_Mesh)
			{
				mesh = (T_Mesh) scene;

				// no bounding box hit? -> next object
				if (!Octree.boundingBoxHit(mesh, eye, ray)) continue;

				float t;
				float[] n;
				float a, rayVn, pen;
				float[] p1, p2, p3;
				float[] ai = new float[3];

				for (int i : object.getValue())
				{
					p1 = mesh.vertices[mesh.triangles[i][0]];
					p2 = mesh.vertices[mesh.triangles[i][1]];
					p3 = mesh.vertices[mesh.triangles[i][2]];

					// fetch precalculated face areas and face normals
					a = mesh.triangleAreas[i];
					n = mesh.triangleNormals[i];

					rayVn = ray[0] * n[0] + ray[1] * n[1] + ray[2] * n[2];

					// backface? => next triangle
					if (rayVn >= 0) continue;

					// no intersection point? => next triangle
					if (Math.abs(rayVn) < EPSILON) continue;

					pen = (p1[0] - eye[0]) * n[0] + (p1[1] - eye[1]) * n[1]
							+ (p1[2] - eye[2]) * n[2];

					// calculate intersection point with plane along the ray
					t = pen / rayVn;

					// already a closer intersection point? => next triangle
					if (t >= minT) continue;

					// the intersection point with the plane
					ip[0] = eye[0] + t * ray[0];
					ip[1] = eye[1] + t * ray[1];
					ip[2] = eye[2] + t * ray[2];
					
					// ### CHANGED ###
					if (!Octree.isPointInsideBox(ip, this.min, this.max)) continue;

					// no intersection point with the triangle? => next
					// triangle
					if (!triangleTest(ip, p1, p2, p3, a, ai)) continue;

					// from here: t < minT and triangle intersection
					// I'm the winner until now!
					minT = t;
					minObjectsKey = scene;
					minIndex = i;

					// prepare everything for shading alternatives
					// the intersection point
					minIP[0] = ip[0];
					minIP[1] = ip[1];
					minIP[2] = ip[2];

					switch (mesh.fgp)
					{
					case 'g':
					case 'G':
						// remember barycentric coordinates bu, bv, bw for shading
						bu = ai[0] / a;
						bv = ai[1] / a;
						bw = ai[2] / a;
						break;
					case 'p':
					case 'P':
						// the normal is barycentrically interpolated between
						// the three vertices
						bu = ai[0] / a;
						bv = ai[1] / a;
						bw = ai[2] / a;
						float nTemp[] = new float[3];
						nTemp[0] = mesh.vertexNormals[mesh.triangles[minIndex][2]][0] * bu
								+ mesh.vertexNormals[mesh.triangles[minIndex][0]][0] * bv
								+ mesh.vertexNormals[mesh.triangles[minIndex][1]][0] * bw;
						nTemp[1] = mesh.vertexNormals[mesh.triangles[minIndex][2]][1] * bu
								+ mesh.vertexNormals[mesh.triangles[minIndex][0]][1] * bv
								+ mesh.vertexNormals[mesh.triangles[minIndex][1]][1] * bw;
						nTemp[2] = mesh.vertexNormals[mesh.triangles[minIndex][2]][2] * bu
								+ mesh.vertexNormals[mesh.triangles[minIndex][0]][2] * bv
								+ mesh.vertexNormals[mesh.triangles[minIndex][1]][2] * bw;
						normalize(nTemp);
						minN = nTemp;

						// the material is barycentrically interpolated between
						// the three vertex materials
						int matIndex0 = mesh.verticesMat[mesh.triangles[minIndex][0]];
						int matIndex1 = mesh.verticesMat[mesh.triangles[minIndex][1]];
						int matIndex2 = mesh.verticesMat[mesh.triangles[minIndex][2]];
						float materialTemp[] = new float[9];
						int materialNTemp;
						for (int k = 0; k < 9; ++k)
						{
							materialTemp[k] = mesh.materials[matIndex2][k] * bu + mesh.materials[matIndex0][k] * bv
									+ mesh.materials[matIndex1][k] * bw;
						}
						minMaterial = materialTemp;
						materialNTemp = (int) (mesh.materialsN[matIndex2] * bu + mesh.materialsN[matIndex0] * bv
								+ mesh.materialsN[matIndex1] * bw);
						minMaterialN = materialNTemp;
					}
				}
			}
		}

		// no intersection point found => return with no result
		if (minObjectsKey == null) return null;

		// light vector at the intersection point
		l[0] = this.centerOfPointLight[0] - minIP[0];
		l[1] = this.centerOfPointLight[1] - minIP[1];
		l[2] = this.centerOfPointLight[2] - minIP[2];
		normalize(l);

		// decide which shading model will be applied
		// implicit: only phong shading available => shade=illuminate
		if (minObjectsKey instanceof I_Sphere)
		{
			return phongIlluminate(minMaterial, minMaterialN, l, minN, v);
		} // triangle mesh: flat, gouraud or phong shading according to file data
		else if (minObjectsKey.getHeader() == "TRIANGLE_MESH")
		{
			mesh = (T_Mesh) minObjectsKey;
			switch (mesh.fgp)
			{
			case 'f':
			case 'F':
				// // illumination can be calculated here
				// // this is a variant between flat und phong shading
				// return phongIlluminate(minMaterial, minMaterialN, l, minN, v, ambientLightColor, specularLight);

				// lookup triangle color of triangle hit
				return new Color(mesh.triangleColors[minIndex][0], mesh.triangleColors[minIndex][1],
						mesh.triangleColors[minIndex][2]);
			case 'g':
			case 'G':
				// the color is barycentrically interpolated between the three
				// vertex colors
				float colorf[] = new float[3];
				colorf[0] = Math.min(mesh.vertexColors[mesh.triangles[minIndex][2]][0] * bu
						+ mesh.vertexColors[mesh.triangles[minIndex][0]][0] * bv
						+ mesh.vertexColors[mesh.triangles[minIndex][1]][0] * bw, 1);
				colorf[1] = Math.min(mesh.vertexColors[mesh.triangles[minIndex][2]][1] * bu
						+ mesh.vertexColors[mesh.triangles[minIndex][0]][1] * bv
						+ mesh.vertexColors[mesh.triangles[minIndex][1]][1] * bw, 1);
				colorf[2] = Math.min(mesh.vertexColors[mesh.triangles[minIndex][2]][2] * bu
						+ mesh.vertexColors[mesh.triangles[minIndex][0]][2] * bv
						+ mesh.vertexColors[mesh.triangles[minIndex][1]][2] * bw, 1);
				return new Color(colorf[0], colorf[1], colorf[2]);
			case 'p':
			case 'P':
				// calculate the color per per pixel phong lightning
				return phongIlluminate(minMaterial, minMaterialN, l, minN, v);
			}
		}
		return null;
	}
	
	private static boolean boundingBoxHit(RT_Object object, float[] eye, float[] ray)
	{
		float t;
		float intersectionPoint[] = new float[3];

		// front and back
		if (Math.abs(ray[2]) > EPSILON)
		{
			// n = [0, 0, 1]
			// front xy
			t = (object.max[2] - eye[2]) / ray[2];

			intersectionPoint[0] = eye[0] + t * ray[0];
			intersectionPoint[1] = eye[1] + t * ray[1];

			if ((intersectionPoint[0] > object.min[0] && intersectionPoint[0] < object.max[0]) &&
				(intersectionPoint[1] > object.min[1] && intersectionPoint[1] < object.max[1])) return true;

			// n = [0, 0, -1]
			// back xy
			t = (object.min[2] - eye[2]) / ray[2];

			intersectionPoint[0] = eye[0] + t * ray[0];
			intersectionPoint[1] = eye[1] + t * ray[1];

			if ((intersectionPoint[0] > object.min[0] && intersectionPoint[0] < object.max[0]) &&
				(intersectionPoint[1] > object.min[1] && intersectionPoint[1] < object.max[1])) return true;
		}

		// left and right
		if (Math.abs(ray[0]) > EPSILON)
		{
			// n = [-1, 0, 0]
			// left xy
			t = (object.min[0] - eye[0]) / ray[0];

			intersectionPoint[1] = eye[1] + t * ray[1];
			intersectionPoint[2] = eye[2] + t * ray[2];

			if ((intersectionPoint[1] > object.min[1] && intersectionPoint[1] < object.max[1]) &&
				(intersectionPoint[2] > object.min[2] && intersectionPoint[2] < object.max[2])) return true;

			// n = [1, 0, 0]
			// right xy
			t = (object.max[0] - eye[0]) / ray[0];

			intersectionPoint[1] = eye[1] + t * ray[1];
			intersectionPoint[2] = eye[2] + t * ray[2];

			if ((intersectionPoint[1] > object.min[1] && intersectionPoint[1] < object.max[1]) &&
				(intersectionPoint[2] > object.min[2] && intersectionPoint[2] < object.max[2])) return true;
		}
		
		// top and bottom
		if (Math.abs(ray[1]) > EPSILON)
		{
			// n = [0, 1, 0]
			// top xy
			t = (object.max[1] - eye[1]) / ray[1];

			intersectionPoint[0] = eye[0] + t * ray[0];
			intersectionPoint[2] = eye[2] + t * ray[2];

			if ((intersectionPoint[0] > object.min[0] && intersectionPoint[0] < object.min[0]) &&
				(intersectionPoint[2] > object.min[2] && intersectionPoint[2] < object.max[2])) return true;

			// n = [0, -1, 0]
			// bottom xy
			t = (object.min[1] - eye[1]) / ray[1];

			intersectionPoint[0] = eye[0] + t * ray[0];
			intersectionPoint[2] = eye[2] + t * ray[2];

			if ((intersectionPoint[0] > object.min[0] && intersectionPoint[0] < object.min[0]) &&
				(intersectionPoint[2] > object.min[2] && intersectionPoint[2] < object.max[2])) return true;
		}
		return false;
	}
	
	private Color phongIlluminate(float[] material, float materialN, float[] l, float[] n, float[] v)
	{
		float ir = 0, ig = 0, ib = 0; // reflected intensity, rgb channels
		float[] reflectionVector = new float[3];
		float scalar_l_n, scalar_r_v;
		scalar_l_n = l[0] * n[0] + l[1] * n[1] + l[2] * n[2];

		// ambient component, ambientLightColor*ra
		ir += this.ambientLightColor[0] * material[0];
		ig += this.ambientLightColor[1] * material[1];
		ib += this.ambientLightColor[2] * material[2];

		// diffuse component, specularLight*rd*<l,n>
		if (scalar_l_n > 0)
		{
			ir += this.specularLight[0] * material[3] * scalar_l_n;
			ig += this.specularLight[1] * material[4] * scalar_l_n;
			ib += this.specularLight[2] * material[5] * scalar_l_n;
			
			// reflection vector r=2*<l,n>*n-l
			reflectionVector[0] = 2f * scalar_l_n * n[0] - l[0];
			reflectionVector[1] = 2f * scalar_l_n * n[1] - l[1];
			reflectionVector[2] = 2f * scalar_l_n * n[2] - l[2];
			normalize(reflectionVector);

			// <r,v>
			scalar_r_v = reflectionVector[0] * v[0] + reflectionVector[1] * v[1] + reflectionVector[2] * v[2];

			// specular component, specularLight*rs*<r,v>^n
			if (scalar_r_v > 0)
			{
				float pow = (float)Math.pow(scalar_r_v, materialN);
				ir += this.specularLight[0] * material[6] * pow;
				ig += this.specularLight[1] * material[7] * pow;
				ib += this.specularLight[2] * material[8] * pow;
			}
		}

		ir = Math.max(Math.min(ir, 1), 0);
		ig = Math.max(Math.min(ig, 1), 0);
		ib = Math.max(Math.min(ib, 1), 0);
		return new Color(ir, ig, ib);
	}
	
	private float calculateIntersectionPointRayFactor(float[] eye, float[] ray)
	{
		float minFactor = Float.POSITIVE_INFINITY;
		float[] intersectionPoint = new float[3];
		if (Octree.isIntersecting(eye, ray, Octree.VOXEL_RIGHT_FACE_NORMAL, this.max, intersectionPoint) &&
			Octree.isPointInsideBox(intersectionPoint, this.min, this.max))
		{
			float factor = Octree.calculateIntersectionPointFactor(eye, ray, intersectionPoint);
			if (0 < factor && minFactor > factor) minFactor = factor;
		}
		if (Octree.isIntersecting(eye, ray, Octree.VOXEL_LEFT_FACE_NORMAL, this.min, intersectionPoint) &&
			Octree.isPointInsideBox(intersectionPoint, this.min, this.max))
		{
			float factor = Octree.calculateIntersectionPointFactor(eye, ray, intersectionPoint);
			if (0 < factor && minFactor > factor) minFactor = factor;
		}
		if (Octree.isIntersecting(eye, ray, Octree.VOXEL_TOP_FACE_NORMAL, this.max, intersectionPoint) &&
			Octree.isPointInsideBox(intersectionPoint, this.min, this.max))
		{
			float factor = Octree.calculateIntersectionPointFactor(eye, ray, intersectionPoint);
			if (0 < factor && minFactor > factor) minFactor = factor;
		}
		if (Octree.isIntersecting(eye, ray, Octree.VOXEL_BOTTOM_FACE_NORMAL, this.min, intersectionPoint) &&
			Octree.isPointInsideBox(intersectionPoint, this.min, this.max))
		{
			float factor = Octree.calculateIntersectionPointFactor(eye, ray, intersectionPoint);
			if (0 < factor && minFactor > factor) minFactor = factor;
		}
		if (Octree.isIntersecting(eye, ray, Octree.VOXEL_FRONT_FACE_NORMAL, this.max, intersectionPoint) &&
			Octree.isPointInsideBox(intersectionPoint, this.min, this.max))
		{
			float factor = Octree.calculateIntersectionPointFactor(eye, ray, intersectionPoint);
			if (0 < factor && minFactor > factor) minFactor = factor;
		}
		if (Octree.isIntersecting(eye, ray, Octree.VOXEL_BACK_FACE_NORMAL, this.min, intersectionPoint) &&
			Octree.isPointInsideBox(intersectionPoint, this.min, this.max))
		{
			float factor = Octree.calculateIntersectionPointFactor(eye, ray, intersectionPoint);
			if (0 < factor && minFactor > factor) minFactor = factor;
		}

		return minFactor;
	}
}