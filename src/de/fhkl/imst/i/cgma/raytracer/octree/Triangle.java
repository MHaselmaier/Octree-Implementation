package de.fhkl.imst.i.cgma.raytracer.octree;

import de.fhkl.imst.i.cgma.raytracer.file.T_Mesh;

final class Triangle
{
	protected final float[] firstPoint;
	protected final float[] secondPoint;
	protected final float[] thirdPoint;
	protected final float[] normal;
	protected final float area;
	
	protected Triangle(T_Mesh mesh, int triangleIndex)
	{
		this.firstPoint = mesh.vertices[mesh.triangles[triangleIndex][0]];
		this.secondPoint = mesh.vertices[mesh.triangles[triangleIndex][1]];
		this.thirdPoint = mesh.vertices[mesh.triangles[triangleIndex][2]];
		this.normal = mesh.triangleNormals[triangleIndex];
		this.area = mesh.triangleAreas[triangleIndex];
	}
}