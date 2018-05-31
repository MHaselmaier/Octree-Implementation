package de.fhkl.imst.i.cgma.raytracer.octree;

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
	private static final int MIN_TRIANGLES = 10;
	
	private static final float[] VOXEL_RIGHT_FACE_NORMAL = {1, 0, 0};
	private static final float[] VOXEL_LEFT_FACE_NORMAL = {-1, 0, 0};
	private static final float[] VOXEL_TOP_FACE_NORMAL = {0, 1, 0};
	private static final float[] VOXEL_BOTTOM_FACE_NORMAL = {0, -1, 0};
	private static final float[] VOXEL_FRONT_FACE_NORMAL = {0, 0, 1};
	private static final float[] VOXEL_BACK_FACE_NORMAL = {0, 0, -1};
	
	private float[] min;
	private float[] max;
	
	public int depth = 0;
	private Octree parent;
	private Octree[] children;
	private Map<RT_Object, Integer[]> objects;
	
	private Octree(Map<RT_Object, Integer[]> objects, float[] min, float[] max)
	{
		this.min = min;
		this.max = max;		
		this.parent = null;
		this.objects = objects;
		generateChildren();
	}
	
	private Octree(Octree parent, float[] min, float[] max)
	{
		this.min = min;
		this.max = max;
		this.parent = parent;
		this.objects = findContainedObjects();
		generateChildren();
	}
	
	private void goDeeper()
	{
		if (null != this.parent)
		{
			this.parent.goDeeper();
		}
		++depth;		
	}
	
	private Map<RT_Object, Integer[]> findContainedObjects()
	{
		Map<RT_Object, Integer[]> containedObjects = new HashMap<>();
		
		for (Map.Entry<RT_Object, Integer[]> object: this.parent.objects.entrySet())
		{
			RT_Object scene = object.getKey();
			
			if (!isBoundingBoxInsideVoxel(scene))
			{
				continue;
			}
			
			if (scene instanceof I_Sphere)
			{
				// Nur Bounding Box Check ...
				// Ist ja eh implizit und besteht nicht aus Dreiecken.
				containedObjects.put(scene, new Integer[1]);
			}
			else if (scene instanceof T_Mesh)
			{
				T_Mesh mesh = (T_Mesh)scene;
				Set<Integer> indices = new TreeSet<>();
				// Teste für jedes Dreieck des Objektes, 
				// ob es innerhalb des Voxels liegt.
				for (int i = 0; mesh.triangles.length > i; ++i)
				{
					if (isTriangleInsideVoxel(mesh.vertices[mesh.triangles[i][0]],
											  mesh.vertices[mesh.triangles[i][1]],
											  mesh.vertices[mesh.triangles[i][2]]))
					{
						indices.add(i);
						continue;
					}					
				}
				if (indices.size() > 0)
				{
					containedObjects.put(mesh, indices.toArray(new Integer[0]));
				}
			}
		}
		return containedObjects;
	}
	
	private boolean isBoundingBoxInsideVoxel(RT_Object scene)
	{
		return (isPointInsideVoxel(scene.min[0], scene.min[1], scene.min[2]) ||
				isPointInsideVoxel(scene.min[0], scene.min[1], scene.max[2]) ||
				isPointInsideVoxel(scene.min[0], scene.max[1], scene.min[2]) ||
				isPointInsideVoxel(scene.min[0], scene.max[1], scene.max[2]) ||
				isPointInsideVoxel(scene.max[0], scene.min[1], scene.min[2]) ||
				isPointInsideVoxel(scene.max[0], scene.min[1], scene.max[2]) ||
				isPointInsideVoxel(scene.max[0], scene.max[1], scene.min[2]) ||
				isPointInsideVoxel(scene.max[0], scene.max[1], scene.max[2]));
	}
	
	private boolean isTriangleInsideVoxel(float[] p1, float[] p2, float[] p3)
	{
		return (isPointInsideVoxel(p1) ||
				isPointInsideVoxel(p2) ||
				isPointInsideVoxel(p3) ||
				isTriangleEdgeIntersectingVoxel(p1, p2, p3));
	}
	
	private boolean isPointInsideVoxel(float[] p)
	{
		return isPointInsideVoxel(p[0], p[1], p[2]);
	}
	
	private boolean isPointInsideVoxel(float x, float y, float z)
	{
		return ((this.min[0] <= x && this.max[0] >= x) &&
				(this.min[1] <= y && this.max[1] >= y) &&
				(this.min[2] <= z && this.max[2] >= z));
	}
	
	private boolean isTriangleEdgeIntersectingVoxel(float[] p1, float[] p2, float[] p3)
	{
		float[] p1p2 = {p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]};
		float[] p1p3 = {p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2]};
		float[] p2p3 = {p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2]};
		
		return (isTriangleEdgeIntersectingFace(p1, p2, p1p2, p1p3, p2p3, VOXEL_RIGHT_FACE_NORMAL, this.max) ||
				isTriangleEdgeIntersectingFace(p1, p2, p1p2, p1p3, p2p3, VOXEL_LEFT_FACE_NORMAL, this.min) ||
				isTriangleEdgeIntersectingFace(p1, p2, p1p2, p1p3, p2p3, VOXEL_TOP_FACE_NORMAL, this.max) ||
				isTriangleEdgeIntersectingFace(p1, p2, p1p2, p1p3, p2p3, VOXEL_BOTTOM_FACE_NORMAL, this.min) ||
				isTriangleEdgeIntersectingFace(p1, p2, p1p2, p1p3, p2p3, VOXEL_FRONT_FACE_NORMAL, this.max) ||
				isTriangleEdgeIntersectingFace(p1, p2, p1p2, p1p3, p2p3, VOXEL_BACK_FACE_NORMAL, this.min));
	}
	
	private boolean isTriangleEdgeIntersectingFace(float[] p1, float[] p2, float[] p1p2, float[] p1p3, float[] p2p3, float[] n, float[] facePoint)
	{
		float[] ip = new float[3];
		return (((intersects(p1, p1p2, n, facePoint, ip) && (1 >= (ip[0] - p1[0]) / p1p2[0])) ||
				 (intersects(p1, p1p3, n, facePoint, ip) && (1 >= (ip[0] - p1[0]) / p1p3[0])) ||
				 (intersects(p2, p2p3, n, facePoint, ip) && (1 >= (ip[0] - p2[0]) / p2p3[0]))) &&
				isPointInsideVoxel(ip));
	}
	
	private boolean intersects(float[] e, float[] v, float[] n, float[] p, float[] ip)
	{
		float vn = v[0] * n[0] + v[1] * n[1] + v[2] * n[2];
		if (Math.abs(vn) < 1E-7)
			return false;
		
		float pen = (p[0] - e[0]) * n[0] + (p[1] - e[1]) * n[1] + (p[2] - e[2]) * n[2];
		float t = pen / vn;

		ip[0] = e[0] + t * v[0];
		ip[1] = e[1] + t * v[1];
		ip[2] = e[2] + t * v[2];
		
		return true;
	}
	
	private void generateChildren()
	{
		if (countContainedTriangles() <= Octree.MIN_TRIANGLES)
		{
			return;
		}
		
		float[] delta = new float[3];
		delta[0] = (this.max[0] - this.min[0]) / 2f;
		delta[1] = (this.max[1] - this.min[1]) / 2f;
		delta[2] = (this.max[2] - this.min[2]) / 2f;
		
		goDeeper();
		
		this.children = new Octree[8];
		float[] newMin;
		float[] newMax;
		
		// Bottom-Back-Left
		newMin = new float[] {this.min[0], this.min[1], this.min[2]};
		newMax = new float[] {this.min[0] + delta[0], this.min[1] + delta[1], this.min[2] + delta[2]};
		this.children[0] = new Octree(this, newMin, newMax);
		
		// Bottom-Back-Right
		newMin = new float[] {this.min[0] + delta[0], this.min[1], this.min[2]};
		newMax = new float[] {this.max[0], this.min[1] + delta[1], this.min[2] + delta[2]};
		this.children[1] = new Octree(this, newMin, newMax);
		
		// Bottom-Front-Right
		newMin = new float[] {this.min[0] + delta[0], this.min[1], this.min[2] + delta[2]};
		newMax = new float[] {this.max[0], this.min[1] + delta[1], this.max[2]};
		this.children[2] = new Octree(this, newMin, newMax);
		
		// Bottom-Front-Left
		newMin = new float[] {this.min[0], this.min[1], this.min[2] + delta[2]};
		newMax = new float[] {this.min[0] + delta[0], this.min[1] + delta[1], this.max[2]};
		this.children[3] = new Octree(this, newMin, newMax);
		
		// Top-Back-Left
		newMin = new float[] {this.min[0], this.min[1] + delta[1], this.min[2]};
		newMax = new float[] {this.min[0] + delta[0], this.max[1], this.min[2] + delta[2]};
		this.children[4] = new Octree(this, newMin, newMax);
		
		// Top-Back-Right
		newMin = new float[] {this.min[0] + delta[0], this.min[1] + delta[1], this.min[2]};
		newMax = new float[] {this.max[0], this.max[1], this.min[2] + delta[2]};
		this.children[5] = new Octree(this, newMin, newMax);
		
		// Top-Front-Right
		newMin = new float[] {this.min[0] + delta[0], this.min[1] + delta[1], this.min[2] + delta[2]};
		newMax = new float[] { this.max[0], this.max[1], this.max[2]};
		this.children[6] = new Octree(this, newMin, newMax);
		
		// Top-Front-Left
		newMin = new float[] {this.min[0], this.min[1] + delta[1], this.min[2] + delta[2]};
		newMax = new float[] {this.min[0] + delta[0], this.max[1], this.max[2]};
		this.children[7] = new Octree(this, newMin, newMax);
	}
	
	private int countContainedTriangles()
	{
		int containedTriangles = 0;
		for (Integer[] indices: this.objects.values())
		{
			containedTriangles += indices.length;
		}
		return containedTriangles;
	}
	
	public static Octree generate(Vector<RT_Object> objects)
	{
		float[] min = new float[3];
		float[] max = new float[3];
		calculateBoundingBox(objects, min, max);
		Map<RT_Object, Integer[]> containedObjects = generateObjectTriangleIndexMap(objects);
		return new Octree(containedObjects, min, max);
	}
	
	private static void calculateBoundingBox(Vector<RT_Object> objects, float[] min, float[] max)
	{
		min[0] = Float.POSITIVE_INFINITY;
		min[1] = Float.POSITIVE_INFINITY;
		min[2] = Float.POSITIVE_INFINITY;
		max[0] = Float.NEGATIVE_INFINITY;
		max[1] = Float.NEGATIVE_INFINITY;
		max[2] = Float.NEGATIVE_INFINITY;
		
		for (RT_Object object: objects)
		{
			if (min[0] > object.min[0]) min[0] = object.min[0];
			if (min[1] > object.min[1]) min[1] = object.min[1];
			if (min[2] > object.min[2]) min[2] = object.min[2];
			
			if (max[0] < object.max[0]) max[0] = object.max[0];
			if (max[1] < object.max[1]) max[1] = object.max[1];
			if (max[2] < object.max[2]) max[2] = object.max[2];
		}
	}
	
	private static Map<RT_Object, Integer[]> generateObjectTriangleIndexMap(Vector<RT_Object> objects)
	{
		Map<RT_Object, Integer[]> containedObjects = new HashMap<>(); 
		for (RT_Object object: objects)
		{
			Integer[] indices = new Integer[1];
			if (object instanceof T_Mesh)
			{
				indices = new Integer[((T_Mesh)object).triangles.length];
			}
			for (int i = 0; indices.length > i; ++i)
			{
				indices[i] = i;
			}
			containedObjects.put(object, indices);
		}
		return containedObjects;
	}
}