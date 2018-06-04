package de.fhkl.imst.i.cgma.raytracer.octree;

import java.awt.Color;
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
	// one hardcoded point light as a minimal solution :-(
	private float[] Ia = { 0.25f, 0.25f, 0.25f }; // ambient light color
	private float[] Ids = { 1.0f, 1.0f, 1.0f }; // diffuse and specular light
	// color
	private float[] ICenter = { 4.0f, 4.0f, 2.0f }; // center of point light
	
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
	
	public Color traceRay(float[] entryPoint, float[] ray)
	{
		if (null == this.children) return traceRayInsideVoxel(entryPoint, ray);
		
		Octree next;
		Color color;
		while (null != entryPoint)
		{
			next = findNextVoxel(entryPoint, ray);
			color = next.traceRay(entryPoint, ray);
			if (null != color)
			{
				return color;
			}
			else
			{
				entryPoint = findExitPoint(next, entryPoint, ray);
			}
		}
		
		return null;
	}
	
	private float[] findExitPoint(Octree next, float[] entryPoint, float[] ray)
	{
		
	}
	
	private Color traceRayInsideVoxel(float[] entryPoint, float[] ray)
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
		v[0] = entryPoint[0] - ray[0];
		v[1] = entryPoint[1] - ray[1];
		v[2] = entryPoint[2] - ray[2];
		normalize(v);

		RT_Object scene;
		I_Sphere sphere;
		T_Mesh mesh;
		
		for (Map.Entry<RT_Object, Integer[]> object: this.objects.entrySet())
		{
			scene = object.getKey();
			if (scene instanceof I_Sphere)
			{
				sphere = (I_Sphere) scene;

				// no bounding box hit? -> next object
				if (!bboxHit(sphere, entryPoint, ray))
					continue;

				// ray intersection uses quadratic equation
				float a, b, c, d;
				a = ray[0] * ray[0] + ray[1] * ray[1] + ray[2] * ray[2];
				b = 2f * (ray[0] * (entryPoint[0] - sphere.center[0]) +
						  ray[1] * (entryPoint[1] - sphere.center[1]) +
						  ray[2] * (entryPoint[2] - sphere.center[2]));
				c = ((entryPoint[0] - sphere.center[0]) * (entryPoint[0] - sphere.center[0]) +
					 (entryPoint[1] - sphere.center[1]) * (entryPoint[1] - sphere.center[1]) +
					 (entryPoint[2] - sphere.center[2]) * (entryPoint[2] - sphere.center[2])) -
					sphere.radius * sphere.radius;

				// positive discriminant determines intersection
				d = b * b - 4f * a * c;
				// no intersection point? => next object
				if (d <= 0)	continue;

				// from here: intersection takes place!

				// calculate first intersection point with sphere along the ray
				float t = (float)(-b - Math.sqrt(d)) / (2 * a);

				// already a closer intersection point? => next object
				if (t >= minT) continue;
				
				ip[0] = entryPoint[0] + t * ray[0];
				ip[1] = entryPoint[1] + t * ray[1];
				ip[2] = entryPoint[2] + t * ray[2];
				if (!isPointInsideVoxel(ip)) continue;

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
				mesh = (T_Mesh)scene;
				
				// no bounding box hit? -> next object
				if (!bboxHit(mesh, entryPoint, ray))
					continue;
				
				float t;
				float[] n;
				float a, rayVn, pen;
				float[] p1, p2, p3;
				float[] ai = new float[3];
				
				for (int i: object.getValue())
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
					if (Math.abs(rayVn) < 1E-7)	continue;

					pen = (p1[0] - entryPoint[0]) * n[0] + (p1[1] - entryPoint[1]) * n[1] + (p1[2] - entryPoint[2]) * n[2];

					// calculate intersection point with plane along the ray
					t = pen / rayVn;

					// already a closer intersection point? => next triangle
					if (t >= minT) continue;

					// the intersection point with the plane
					ip[0] = entryPoint[0] + t * ray[0];
					ip[1] = entryPoint[1] + t * ray[1];
					ip[2] = entryPoint[2] + t * ray[2];
					if (!isPointInsideVoxel(ip)) continue;

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
						nTemp[0] = mesh.vertexNormals[mesh.triangles[minIndex][2]][0] * bu +
								   mesh.vertexNormals[mesh.triangles[minIndex][0]][0] * bv +
								   mesh.vertexNormals[mesh.triangles[minIndex][1]][0] * bw;
						nTemp[1] = mesh.vertexNormals[mesh.triangles[minIndex][2]][1] * bu +
								   mesh.vertexNormals[mesh.triangles[minIndex][0]][1] * bv +
								   mesh.vertexNormals[mesh.triangles[minIndex][1]][1] * bw;
						nTemp[2] = mesh.vertexNormals[mesh.triangles[minIndex][2]][2] * bu +
								   mesh.vertexNormals[mesh.triangles[minIndex][0]][2] * bv +
								   mesh.vertexNormals[mesh.triangles[minIndex][1]][2] * bw;
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
							materialTemp[k] = mesh.materials[matIndex2][k] * bu +
											  mesh.materials[matIndex0][k] * bv +
											  mesh.materials[matIndex1][k] * bw;
						}
						minMaterial = materialTemp;
						materialNTemp = (int)(mesh.materialsN[matIndex2] * bu +
											  mesh.materialsN[matIndex0] * bv +
											  mesh.materialsN[matIndex1] * bw);
						minMaterialN = materialNTemp;
					}
				}
			}
		}
		
		// no intersection point found => return with no result
		if (minObjectsKey == null) return null;

		// light vector at the intersection point
		l[0] = ICenter[0] - minIP[0];
		l[1] = ICenter[1] - minIP[1];
		l[2] = ICenter[2] - minIP[2];
		normalize(l);

		// decide which shading model will be applied

		// implicit: only phong shading available => shade=illuminate
		if (minObjectsKey instanceof I_Sphere)
		{
			return phongIlluminate(minMaterial, minMaterialN, l, minN, v);
		}
		// triangle mesh: flat, gouraud or phong shading according to file data
		else if (minObjectsKey.getHeader() == "TRIANGLE_MESH")
		{
			mesh = (T_Mesh)minObjectsKey;
			switch (mesh.fgp)
			{
			case 'f':
			case 'F':
				// // illumination can be calculated here
				// // this is a variant between flat und phong shading
				// return phongIlluminate(minMaterial, minMaterialN, l, minN, v, Ia, Ids);

				// lookup triangle color of triangle hit
				return new Color(mesh.triangleColors[minIndex][0],
						         mesh.triangleColors[minIndex][1],
						         mesh.triangleColors[minIndex][2]);
			case 'g':
			case 'G':
				// the color is barycentrically interpolated between the three
				// vertex colors
				float colorf[] = new float[3];
				colorf[0] = Math.min(mesh.vertexColors[mesh.triangles[minIndex][2]][0] * bu	+
						             mesh.vertexColors[mesh.triangles[minIndex][0]][0] * bv +
						             mesh.vertexColors[mesh.triangles[minIndex][1]][0] * bw, 1);
				colorf[1] = Math.min(mesh.vertexColors[mesh.triangles[minIndex][2]][1] * bu	+
						             mesh.vertexColors[mesh.triangles[minIndex][0]][1] * bv	+
						             mesh.vertexColors[mesh.triangles[minIndex][1]][1] * bw, 1);
				colorf[2] = Math.min(mesh.vertexColors[mesh.triangles[minIndex][2]][2] * bu	+
						             mesh.vertexColors[mesh.triangles[minIndex][0]][2] * bv +
						             mesh.vertexColors[mesh.triangles[minIndex][1]][2] * bw, 1);
				return new Color(colorf[0], colorf[1], colorf[2]);
			case 'p':
			case 'P':
				// calculate the color per per pixel phong lightning
				return phongIlluminate(minMaterial, minMaterialN, l, minN, v);
			}
		}
		return null;
	}
	
	private boolean bboxHit(RT_Object object, float[] entryPoint, float[] ray)
	{
		float t;
		float ip[] = new float[3];

		// front and back
		if (Math.abs(ray[2]) > 1E-5)
		{
			// n = [0, 0, 1]
			// front xy
			t = (object.max[2] - entryPoint[2]) / ray[2];

			ip[0] = entryPoint[0] + t * ray[0];
			ip[1] = entryPoint[1] + t * ray[1];

			if (ip[0] > object.min[0] && ip[0] < object.max[0] && ip[1] > object.min[1] && ip[1] < object.max[1])
			{
				return true;
			}

			// n = [0, 0, -1]
			// back xy
			t = (entryPoint[2] - object.min[2]) / ray[2];

			ip[0] = entryPoint[0] + t * ray[0];
			ip[1] = entryPoint[1] + t * ray[1];

			if (ip[0] > object.min[0] && ip[0] < object.max[0] && ip[1] > object.min[1] && ip[1] < object.max[1])
			{
				return true;
			}
		}

		// left and right
		if (Math.abs(ray[0]) > 1E-5)
		{
			// n = [-1, 0, 0]
			// left xy
			t = (entryPoint[0] - object.min[0]) / ray[0];

			ip[1] = entryPoint[1] + t * ray[1];
			ip[2] = entryPoint[2] + t * ray[2];

			if (ip[1] > object.min[1] && ip[1] < object.max[1] && ip[2] > object.min[2] && ip[2] < object.max[2])
			{
				return true;
			}

			// n = [1, 0, 0]
			// right xy
			t = (object.max[0] - entryPoint[0]) / ray[0];

			ip[1] = entryPoint[1] + t * ray[1];
			ip[2] = entryPoint[2] + t * ray[2];

			if (ip[1] > object.min[1] && ip[1] < object.max[1] && ip[2] > object.min[2] && ip[2] < object.max[2])
			{
				return true;
			}
		}
		// top and bottom
		if (Math.abs(ray[1]) > 1E-5)
		{
			// n = [0, 1, 0]
			// top xy
			t = (object.max[1] - entryPoint[1]) / ray[1];

			ip[0] = entryPoint[0] + t * ray[0];
			ip[2] = entryPoint[2] + t * ray[2];

			if (ip[0] > object.min[0] && ip[0] < object.min[0] && ip[2] > object.min[2] && ip[2] < object.max[2])
			{
				return true;
			}

			// n = [0, -1, 0]
			// bottom xy
			t = (entryPoint[1] - object.min[1]) / ray[1];

			ip[0] = entryPoint[0] + t * ray[0];
			ip[2] = entryPoint[2] + t * ray[2];

			if (ip[0] > object.min[0] && ip[0] < object.min[0] && ip[2] > object.min[2] && ip[2] < object.max[2])
			{
				return true;
			}
		}
		return false;
	}
	
	private boolean triangleTest(float[] p, float[] p1, float[] p2, float[] p3, float a, float ai[])
	{
		ai[0] = calculateArea(p1, p2, p);
		ai[1] = calculateArea(p, p2, p3);
		ai[2] = calculateArea(p1, p, p3);

		return (1E-5 > Math.abs(a - (ai[0] + ai[1] + ai[2])));
	}
	
	private Color phongIlluminate(float[] material, float materialN, float[] l, float[] n, float[] v)
	{
		float ir = 0, ig = 0, ib = 0; // reflected intensity, rgb channels
		float[] r = new float[3]; // reflection vector
		float ln, rv; // scalar products <l,n> and <r,v>

		// <l,n>
		ln = l[0] * n[0] + l[1] * n[1] + l[2] * n[2];

		// ambient component, Ia*ra
		ir += this.Ia[0] * material[0];
		ig += this.Ia[1] * material[1];
		ib += this.Ia[2] * material[2];

		// diffuse component, Ids*rd*<l,n>
		if (ln > 0)
		{
			ir += this.Ids[0] * material[3] * ln;
			ig += this.Ids[1] * material[4] * ln;
			ib += this.Ids[2] * material[5] * ln;

			// reflection vector r=2*<l,n>*n-l
			r[0] = 2f * ln * n[0] - l[0];
			r[1] = 2f * ln * n[1] - l[1];
			r[2] = 2f * ln * n[2] - l[2];
			normalize(r);

			// <r,v>
			rv = r[0] * v[0] + r[1] * v[1] + r[2] * v[2];

			// specular component, Ids*rs*<r,v>^n
			if (rv > 0)
			{
				float pow = (float)Math.pow(rv, materialN);
				ir += this.Ids[0] * material[6] * pow;
				ig += this.Ids[1] * material[7] * pow;
				ib += this.Ids[2] * material[8] * pow;
			}
		}

		ir = Math.min(ir, 1);
		ig = Math.min(ig, 1);
		ib = Math.min(ib, 1);
		return new Color(ir, ig, ib);
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
											  mesh.vertices[mesh.triangles[i][2]],
											  mesh.triangleNormals[i],
											  mesh.triangleAreas[i]))
					{
						indices.add(i);
						continue;
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
	
	private boolean isTriangleInsideVoxel(float[] p1, float[] p2, float[] p3, float[] n, float a)
	{
		return (isPointInsideVoxel(p1) ||
				isPointInsideVoxel(p2) ||
				isPointInsideVoxel(p3) ||
				isTriangleEdgeIntersectingVoxel(p1, p2, p3) ||
				isVoxelEdgeIntersectingTriangle(p1, p2, p3, n, a));
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
		return ((isIntersectingInRange(p1, p1p2, n, facePoint, ip) && isPointInsideVoxel(ip)) ||
				(isIntersectingInRange(p1, p1p3, n, facePoint, ip) && isPointInsideVoxel(ip))||
				(isIntersectingInRange(p2, p2p3, n, facePoint, ip) && isPointInsideVoxel(ip)));
	}
	
	private boolean isIntersecting(float[] e, float[] v, float[] n, float[] p, float[] ip)
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
	
	private boolean isIntersectingInRange(float[] e, float[] v, float[] n, float[] p, float[] ip)
	{
		if (isIntersecting(e, v, n, p, ip))
		{
			if (0 < Math.abs(v[0]))
			{
				return (1 >= (ip[0] - e[0]) / v[0]);
			}
			if (0 < Math.abs(v[1]))
			{
				return (1 >= (ip[1] - e[1]) / v[1]);
			}
			if (0 < Math.abs(v[2]))
			{
				return (1 >= (ip[2] - e[2]) / v[2]);
			}
		}
		
		return false;
	}

	private boolean isVoxelEdgeIntersectingTriangle(float[] p1, float[] p2, float[] p3, float[] n, float a)
	{
		float[] bottomBackLeft = {this.min[0], this.min[1], this.min[2]};
		float[] bottomFrontRight = {this.max[0], this.min[1], this.max[2]};
		float[] topBackRight = {this.max[0], this.max[1], this.min[2]};
		float[] topFrontLeft = {this.min[0], this.max[1], this.max[2]};
		
		return (isIntersectingTriangle(bottomBackLeft, new float[] {this.max[0] - this.min[0], 0, 0}, p1, p2, p3, n, a) ||
				isIntersectingTriangle(bottomBackLeft, new float[] {0, this.max[1] - this.min[1], 0}, p1, p2, p3, n, a) ||
				isIntersectingTriangle(bottomBackLeft, new float[] {0, 0, this.max[2] - this.min[2]}, p1, p2, p3, n, a) ||
				isIntersectingTriangle(bottomFrontRight, new float[] {this.min[0] - this.max[0], 0, 0}, p1, p2, p3, n, a) ||
				isIntersectingTriangle(bottomFrontRight, new float[] {0, this.max[1] - this.min[1], 0}, p1, p2, p3, n, a) ||
				isIntersectingTriangle(bottomFrontRight, new float[] {0, 0, this.min[2] - this.max[2]}, p1, p2, p3, n, a) ||
				isIntersectingTriangle(topBackRight, new float[] {this.min[0] - this.max[0], 0, 0}, p1, p2, p3, n, a) ||
				isIntersectingTriangle(topBackRight, new float[] {0, this.min[1] - this.max[1], 0}, p1, p2, p3, n, a) ||
				isIntersectingTriangle(topBackRight, new float[] {0, 0, this.max[2] - this.min[2]}, p1, p2, p3, n, a) ||
				isIntersectingTriangle(topFrontLeft, new float[] {this.max[0] - this.min[0], 0, 0}, p1, p2, p3, n, a) ||
				isIntersectingTriangle(topFrontLeft, new float[] {0, this.min[1] - this.max[1], 0}, p1, p2, p3, n, a) ||
				isIntersectingTriangle(topFrontLeft, new float[] {0, 0, this.min[2] - this.max[2]}, p1, p2, p3, n, a));
	}
	
	private boolean isIntersectingTriangle(float[] e, float[] v, float[] p1, float[] p2, float[] p3, float[] n, float a)
	{
		float[] ip = new float[3];
		return (isIntersectingInRange(e, v, n, p1, ip) && triangleTest(ip, p1, p2, p3, a, new float[3]));
	}
	
	private float calculateArea(float[] p1, float[] p2, float[] p3) {
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
	
	private float normalize(float[] vec) {
		float l = (float)Math.sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
		vec[0] /= l;
		vec[1] /= l;
		vec[2] /= l;
		return l;
	}
	
	private void generateChildren()
	{
		if (Octree.MIN_TRIANGLES >= countContainedTriangles()) return;
		
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