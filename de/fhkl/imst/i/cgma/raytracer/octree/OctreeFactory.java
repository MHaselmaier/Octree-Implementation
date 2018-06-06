package de.fhkl.imst.i.cgma.raytracer.octree;

import de.fhkl.imst.i.cgma.raytracer.file.I_Sphere;
import de.fhkl.imst.i.cgma.raytracer.file.RT_Object;
import de.fhkl.imst.i.cgma.raytracer.file.T_Mesh;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Vector;

public abstract class OctreeFactory
{
    public static Octree buildOctree(Vector<RT_Object> objects)
    {
        float[] min = new float[3];
        float[] max = new float[3];
        OctreeFactory.calculateBoundingBox(objects, min, max);
        return new Octree(OctreeFactory.generateObjectTriangleIndexMap(objects), min, max);
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

}
