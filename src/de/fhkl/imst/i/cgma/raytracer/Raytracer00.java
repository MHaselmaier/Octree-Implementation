package de.fhkl.imst.i.cgma.raytracer;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.Random;
import java.util.Vector;

import de.fhkl.imst.i.cgma.raytracer.file.I_Sphere;
import de.fhkl.imst.i.cgma.raytracer.file.RTFile;
import de.fhkl.imst.i.cgma.raytracer.file.RTFileReader;
import de.fhkl.imst.i.cgma.raytracer.file.RT_Object;
import de.fhkl.imst.i.cgma.raytracer.file.T_Mesh;
import de.fhkl.imst.i.cgma.raytracer.gui.IRayTracerImplementation;
import de.fhkl.imst.i.cgma.raytracer.gui.RayTracerGui;
import de.fhkl.imst.i.cgma.raytracer.octree.Octree;
import java.util.Optional;

public class Raytracer00 implements IRayTracerImplementation {
        // viewing volume with infinite end

        private float fovyDegree;
        private float near;
        private float fovyRadians;

        private float[] Ia = {0.25f, 0.25f, 0.25f}; //ambient
        private float[] Ids = {1.0f, 1.0f, 1.0f};
        private float[] ICenter = {4.0f, 4.0f, 2.0f};

        RayTracerGui gui = new RayTracerGui(this);

        private int resx, resy; // viewport resolution
        private float h, w, aspect; // window height, width and aspect ratio
        
        public boolean isOctreeEnabled = false;

        Vector<RT_Object> objects;

        private Raytracer00() {
                try {
                        gui.addObject(RTFileReader.read(I_Sphere.class, new File("data/ikugel2.dat")));
                        gui.addObject(RTFileReader.read(T_Mesh.class, new File("data/dreiecke2.dat")));
                        gui.addObject(RTFileReader.read(T_Mesh.class, new File("data/kugel1.dat")));
                        gui.addObject(RTFileReader.read(T_Mesh.class, new File("data/kugel2.dat")));
                        gui.addObject(RTFileReader.read(T_Mesh.class, new File("data/kugel3.dat")));

                        objects = gui.getObjects();

                } catch (IOException e) {
                        e.printStackTrace();
                }
        }

        public void setViewParameters(float fovyDegree, float near) {
                // set attributes fovyDegree, fovyRadians, near
                this.fovyDegree = fovyDegree;
                this.fovyRadians = (float) Math.toRadians(fovyDegree);
                this.near = near;

                // set attributes resx, resy, aspect
                resx = gui.getResX();
                resy = gui.getResY();

                aspect = resx / (float) resy;

                // set attributes h, w
                h = (2 * near * (float) Math.tan(this.fovyRadians / 2));
                w = h * aspect;

        }

        @Override
        public void doRayTrace() {
                float x, y, z;			// intersection point in viewing plane
                float rayEx, rayEy, rayEz;	// eye point==ray starting point
                float rayVx, rayVy, rayVz;	// ray vector
                Color color;
          

                // prepare mesh data (normals and areas)
                prepareMeshData();

                Optional<Octree> octree = Optional.empty();
                if(isOctreeEnabled) 
                    octree = Optional.of(new Octree(this.objects));
                
                System.out.println("Octree status:" + isOctreeEnabled);
                
                // hardcoded viewing volume with fovy and near
                setViewParameters(90.0f, 1.0f);
                // set eye point
                rayEx = rayEy = rayEz = 0;

                z = -near;

                //Prepare mesh data for shading
                precalculateMeshDataShading();

                Random rd = new Random();
                // xp, yp: pixel coordinates
                for (int xp = 0; xp < resx; ++xp) {
                        for (int yp = 0; yp < resy; ++yp) {

                                // x, y: view coordinates
                                x = ((xp / (float) (resx - 1)) * w - (w / 2));
                                y = ((resy - 1 - yp) / (float) (resy - 1) * h - (h / 2));

                                // ray vector
                                rayVx = x - rayEx;
                                rayVy = y - rayEy;
                                rayVz = z - rayEz;

                                // get color or null along the ray
                                //color=traceRayAndGetColor...
                                
                                if(octree.isPresent())
                                    color = octree.get().traceRay(new float[]{rayEx, rayEy, rayEz}, new float[]{rayVx, rayVy, rayVz});
                                else
                                    color = traceRayAndGetColor(rayEx, rayEy, rayEz, rayVx, rayVy, rayVz);
                                // if color!=null set pixel with color

                                if (color != null) {
                                        gui.setPixel(xp, yp, color.getRGB());
                                }

                        }
                }
        }

        // returns Color object or null if no intersection was found
        private Color traceRayAndGetColor(float rayEx, float rayEy, float rayEz, float rayVx, float rayVy, float rayVz) {
                // RTFile scene = gui.getFile();

                double minT = Float.MAX_VALUE;
                int minObjectsIndex = -1;
                int minIndex = -1;
                float[] minIP = new float[3];
                float[] minN = new float[3];
                float[] minMaterial = new float[3];
                float minMaterialN = 1;
                float bu = 0, bv = 0, bw = 1;

                float[] v = new float[3];
                float[] l = new float[3];

                //v[0] = rayEx - rayVx;
                //v[1] = rayEy - rayVy;
                //v[2] = rayEz - rayVz;
                v[0] = -rayVx;
                v[1] = -rayVy;
                v[2] = -rayVz;

                normalize(v);

                RTFile scene;
                I_Sphere sphere;
                T_Mesh mesh;

                // loop over all scene objects to find the nearest intersection, that
                // is:
                // object with number minObjectIndex
                // minT is the minimal factor t of the ray equation s(t)=rayE+t*rayV
                // where the nearest intersection takes place
                for (int objectsNumber = 0; objectsNumber < objects.size(); objectsNumber++) {
                        scene = objects.get(objectsNumber);

                        // object is an implicit sphere?
                        if (scene instanceof I_Sphere) {
                                sphere = (I_Sphere) scene;

                                float t;

                                if (!bboxHit(sphere, rayEx, rayEy, rayEz, rayVx, rayVy, rayVz)) {
                                        continue;
                                }

                                // ray intersection uses quadratic equation
                                float a, b, c, d;

                                a = rayVx * rayVx + rayVy * rayVy + rayVz * rayVz;
                                b = (2 * rayVx) * (rayEx - sphere.center[0])
                                        + (2 * rayVy) * (rayEy - sphere.center[1])
                                        + (2 * rayVz) * (rayEz - sphere.center[2]);
                                float cHelper = (rayEx - sphere.center[0]) * (rayEx - sphere.center[0])
                                        + (rayEy - sphere.center[1]) * (rayEy - sphere.center[1])
                                        + (rayEz - sphere.center[2]) * (rayEz - sphere.center[2]);
                                c = cHelper - (sphere.radius * sphere.radius);

                                // positive discriminant determines intersection
                                //d = -42;
                                d = b * b - (4 * a * c);
                                // no intersection point? => next object
                                if (d <= 0) {
                                        continue;
                                }

                                // from here: intersection takes place!
                                // calculate first intersection point with sphere along the
                                // ray
                                t = (float) ((-b - Math.sqrt(d)) / (2 * a));

                                // already a closer intersection point? => next object
                                if (t >= minT) {
                                        continue;
                                }

                                // from here: t < minT
                                // I'm the winner until now!
                                minT = t;
                                minObjectsIndex = objectsNumber;

                                //intersection
                                minIP[0] = rayEx + t * rayVx;
                                minIP[1] = rayEy + t * rayVy;
                                minIP[2] = rayEz + t * rayVz;

                                //normal vector at intersection
                                minN[0] = minIP[0] - sphere.center[0];
                                minN[1] = minIP[1] - sphere.center[1];
                                minN[2] = minIP[2] - sphere.center[2];

                                normalize(minN);

                                //material
                                minMaterial = sphere.material;
                                minMaterialN = sphere.materialN;
                        } else if (scene instanceof T_Mesh) {
                                mesh = (T_Mesh) scene;

                                float t;
                                float[] n;
                                float[] ip = new float[3];

                                //no bounding box hit? -> next object
                                if (!bboxHit(mesh, rayEx, rayEy, rayEz, rayVx, rayVy, rayVz)) {
                                        continue;
                                }

                                float a, rayVn, pen;
                                float[] p1, p2, p3;
                                float[] ai = new float[3];

                                for (int i = 0; i < mesh.triangles.length; i++) {

                                        p1 = mesh.vertices[mesh.triangles[i][0]];
                                        p2 = mesh.vertices[mesh.triangles[i][1]];
                                        p3 = mesh.vertices[mesh.triangles[i][2]];

                                        //a = calculateN(n, p1, p2, p3);
                                        a = mesh.triangleAreas[i];
                                        n = mesh.triangleNormals[i];

                                        rayVn = rayVx * n[0] + rayVy * n[1] + rayVz * n[2];

                                        //backface?
                                        if (rayVn >= 0) {
                                                continue;
                                        }

                                        //no intersection point?
                                        if (Math.abs(rayVn) < 1E-7) {
                                                continue;
                                        }

                                        pen = (p1[0] - rayEx) * n[0]
                                                + (p1[1] - rayEy) * n[1]
                                                + (p1[2] - rayEz) * n[2];

                                        t = pen / rayVn;

                                        if (t >= minT) {
                                                continue;
                                        }

                                        ip[0] = rayEx + t * rayVx;
                                        ip[1] = rayEy + t * rayVy;
                                        ip[2] = rayEz + t * rayVz;

                                        if (!triangleTest(ip, p1, p2, p3, a, ai)) {
                                                continue;
                                        }

                                        minT = t;
                                        minObjectsIndex = objectsNumber;
                                        minIndex = i;

                                        minIP[0] = ip[0];
                                        minIP[1] = ip[1];
                                        minIP[2] = ip[2];

                                        switch (mesh.fgp) {
                                                case 'f':
                                                case 'F':

                                                       /* minN[0] = n[0];
                                                        minN[1] = n[1];
                                                        minN[2] = n[2];

                                                        int matIndex = mesh.verticesMat[mesh.triangles[minIndex][0]];
                                                        minMaterial = mesh.materials[matIndex];
                                                        minMaterialN = mesh.materialsN[matIndex];*/
                                                        break;

                                                case 'g':
                                                case 'G':
                                                        //barycentric coordinates for shading
                                                        bu = ai[0] / a;
                                                        bv = ai[1] / a;
                                                        bw = ai[2] / a;

                                                        break;

                                                case 'p':
                                                case 'P':

                                                        //normal is interpolated between the three vertices
                                                        bu = ai[0] / a;
                                                        bv = ai[1] / a;
                                                        bw = ai[2] / a;

                                                        float nTemp[] = new float[3];
                                                        nTemp[0] = bu * mesh.vertexNormals[mesh.triangles[minIndex][2]][0] 
                                                                + bv * mesh.vertexNormals[mesh.triangles[minIndex][0]][0]
                                                                + bw * mesh.vertexNormals[mesh.triangles[minIndex][1]][0];
                                                        
                                                        nTemp[1] = bu * mesh.vertexNormals[mesh.triangles[minIndex][2]][1] 
                                                                + bv * mesh.vertexNormals[mesh.triangles[minIndex][0]][1]
                                                                + bw * mesh.vertexNormals[mesh.triangles[minIndex][1]][1];
                                                        
                                                        nTemp[2] = bu * mesh.vertexNormals[mesh.triangles[minIndex][2]][2] 
                                                                + bv * mesh.vertexNormals[mesh.triangles[minIndex][0]][2]
                                                                + bw * mesh.vertexNormals[mesh.triangles[minIndex][1]][2];
                                                        
                                                        //normalize(nTemp);

                                                        minN = nTemp;

                                                        // intermediate version
                                                        // the material is not interpolated
                                                        // matIndex =
                                                        // mesh.verticesMat[mesh.triangles[minIndex][0]];
                                                        // minMaterial = mesh.materials[matIndex];
                                                        // minMaterialN = mesh.materialsN[matIndex];
                                                        // the material is barycentrically interpolated between
                                                        // the three vertex materials
                                                        int matIndex0 = mesh.verticesMat[mesh.triangles[minIndex][0]];
                                                        int matIndex1 = mesh.verticesMat[mesh.triangles[minIndex][1]];
                                                        int matIndex2 = mesh.verticesMat[mesh.triangles[minIndex][2]];
                                                        float materialTemp[] = new float[9];
                                                        int materialNTemp;

                                                        for (int k = 0; k < 9; k++) {
                                                                materialTemp[k] = bu * mesh.materials[matIndex2][k] 
                                                                        + bv * mesh.materials[matIndex0][k] 
                                                                        + bw * mesh.materials[matIndex1][k];
                                                        }

                                                        minMaterial = materialTemp;
                                                        materialNTemp = (int)(bu * mesh.materialsN[matIndex2] 
                                                                + bv * mesh.materialsN[matIndex0]
                                                                + bw * mesh.materialsN[matIndex1]);
                                                        minMaterialN = materialNTemp;
                                        }
                                }
                        } else {
                                continue;
                        }
                }

                // no intersection point found => return with no result
                if (minObjectsIndex == -1) {
                        return null;
                }

                // intermediate version
                //Random rd = new Random();
                //return new Color(rd.nextFloat(), rd.nextFloat(), rd.nextFloat());
                //light vector at intersection
                l[0] = ICenter[0] - minIP[0];
                l[1] = ICenter[1] - minIP[1];
                l[2] = ICenter[2] - minIP[2];

                normalize(l);

                if (objects.get(minObjectsIndex) instanceof I_Sphere) {
                        return phongIlluminate(minMaterial, minMaterialN, l, minN, v, Ia, Ids);
                } else if (objects.get(minObjectsIndex).getHeader().equals("TRIANGLE_MESH")) {
                        mesh = ((T_Mesh) objects.get(minObjectsIndex));
                        switch (mesh.fgp) {
                                case 'f':
                                case 'F':

                                        //return phongIlluminate(minMaterial, minMaterialN, l, minN, v, Ia, Ids);
                                        //lookup triangle color of triangle hit
                                        
                                        return new Color(mesh.triangleColors[minIndex][0], mesh.triangleColors[minIndex][1], mesh.triangleColors[minIndex][2]);

                                case 'g':
                                case 'G':
                                        //the color is barycentrically interpolated between the three vertex colors
                                        float colorf[] = new float[3];
                                        colorf[0] = bu * mesh.vertexColors[mesh.triangles[minIndex][2]][0]
                                                + bv * mesh.vertexColors[mesh.triangles[minIndex][0]][0]
                                                + bw * mesh.vertexColors[mesh.triangles[minIndex][1]][0];
                                        
                                        colorf[1] = bu * mesh.vertexColors[mesh.triangles[minIndex][2]][1]
                                                + bv * mesh.vertexColors[mesh.triangles[minIndex][0]][1]
                                                + bw * mesh.vertexColors[mesh.triangles[minIndex][1]][1];
                                        
                                        colorf[2] = bu * mesh.vertexColors[mesh.triangles[minIndex][2]][2]
                                                + bv * mesh.vertexColors[mesh.triangles[minIndex][0]][2]
                                                + bw * mesh.vertexColors[mesh.triangles[minIndex][1]][2];

                                        return new Color(colorf[0] > 1 ? 1 : colorf[0], 
                                                colorf[1] > 1 ? 1 : colorf[1],
                                                colorf[2] > 1 ? 1 : colorf[2]);

                                case 'p':
                                case 'P':
                                        //calculate the color per pixel phong lighting
                                        return phongIlluminate(minMaterial, minMaterialN, l, minN, v, Ia, Ids);
                        }
                }

                return null;
        }

        private Color phongIlluminate(float[] material, float materialN, float[] l, float[] n, float[] v, float[] Ia, float[] Ids) {
                float ir = 0, ig = 0, ib = 0;
                float[] r = new float[3];
                float ln, rv; // scalar products <l,n> and <r,v>

                // <l,n>
                ln = l[0] * n[0] + l[1] * n[1] + l[2] * n[2];

                //Material: 0-2ambient, 3-5 diffuse, 6-8 specular
                // ambient component, Ia*ra
                ir += Ia[0] * material[0];
                ig += Ia[1] * material[1];
                ib += Ia[2] * material[2];

                // diffuse component, Ids*rd*<l,n>
                if (ln > 0) {
                        ir += Ids[0] * material[3] * ln;
                        ig += Ids[1] * material[4] * ln;
                        ib += Ids[2] * material[5] * ln;

                        // reflection vector r=2*<l,n>*n-l
                        r[0] = 2 * ln * n[0] - l[0];
                        r[1] = 2 * ln * n[1] - l[1];
                        r[2] = 2 * ln * n[2] - l[2];

                        normalize(r);

                        // <r,v>
                        rv = r[0] * v[0] + r[1] * v[1] + r[2] * v[2];

                        // specular component, Ids*rs*<r,v>^n
                        if (rv > 0) {
                                float pow = (float) Math.pow(rv, materialN);

                                ir += Ids[0] * material[6] * pow;
                                ir = Math.min(1, ir);
                                ir = Math.max(0, ir);

                                ig += Ids[1] * material[7] * pow;
                                ig = Math.min(1, ig);
                                ig = Math.max(0, ig);

                                ib += Ids[2] * material[8] * pow;
                                ib = Math.min(1, ib);
                                ib = Math.max(0, ib);
                        }
                }

                return new Color(ir, ig, ib);
        }

        // vector normalization
        // CAUTION: vec is an in-/output parameter; the referenced object will be
        // altered!
        private float normalize(float[] vec) {
                float l;

                l = ((float) Math.sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]));
                vec[0] /= l;
                vec[1] /= l;
                vec[2] /= l;

                return l;
        }

        // calculate triangle test
        // is p (the intersection point with the plane through p1, p2 and p3) inside
        // the triangle p1, p2 and p3?
        // the return value answers this question
        // a is an input parameter - the given area of the triangle p1, p2 and p3
        // ai will be computed to be the areas of the sub-triangles to allow to
        // compute barycentric coordinates of the intersection point p
        // ai[0] is associated with bu (p1p2p) across from p3
        // ai[1] is associated with bv (pp2p3) across from p1
        // ai[2] is associated with bw (p1pp3) across form p2
        // CAUTION: ai is an output parameter; the referenced object will be
        // altered!
        private boolean triangleTest(float[] p, float[] p1, float[] p2, float[] p3, float a, float ai[]) {
                float tmp[] = new float[3];

                ai[0] = calculateN(tmp, p1, p2, p);
                ai[1] = calculateN(tmp, p, p2, p3);
                ai[2] = calculateN(tmp, p1, p, p3);

                if (1E-5 > Math.abs((ai[0] + ai[1] + ai[2]) - a)) {
                        return true;
                }

                return false;
        }

        private boolean triangleTest2(float[] p, float[] p1, float[] p2, float[] p3, float a, float ai[]) {
                float tmp[] = new float[3];

                calculateN2(tmp, p1, p2, p);
                calculateN2(tmp, p, p2, p3);
                calculateN2(tmp, p1, p, p3);

                if (1E-5 > Math.abs((ai[0] + ai[1] + ai[2]) - a)) {
                        return true;
                }

                return false;
        }

        // calculate bounding box test
        // decides whether the ray s(t)=rayE+t*rayV intersects the axis aligned
        // bounding box of object -> return value true
        // six plane intersections with rectangle inside tests; if one succeeds
        // bounding box is hit
        private boolean bboxHit(RT_Object object, float rayEx, float rayEy, float rayEz, float rayVx, float rayVy, float rayVz) {
                float t;
                float ip[] = new float[3];

                // front and back
                if (Math.abs(rayVz) > 1E-5) {
                        // front xy
                        t = (object.max[2] - rayEz) / rayVz;
                        ip[0] = rayEx + t * rayVx;
                        ip[1] = rayEy + t * rayVy;

                        if (ip[0] > object.min[0] && ip[0] < object.max[0] && ip[1] > object.min[1] && ip[1] < object.max[1]) {
                                return true;
                        }

                        // back xy
                        t = (object.min[2] - rayEz) / rayVz;
                        ip[0] = rayEx + t * rayVx;
                        ip[1] = rayEy + t * rayVy;

                        if (ip[0] > object.min[0] && ip[0] < object.max[0] && ip[1] > object.min[1] && ip[1] < object.max[1]) {
                                return true;
                        }
                }

                // left and right
                if (Math.abs(rayVx) > 1E-5) {
                        //left yz
                        t = (object.min[0] - rayEx) / rayVx;
                        ip[1] = rayEy + t * rayVy;
                        ip[2] = rayEz + t * rayVz;

                        if (ip[1] < object.max[1] && ip[1] > object.min[1] && ip[2] < object.max[2] && ip[2] > object.min[2]) {
                                return true;
                        }

                        //right yz
                        t = (object.max[0] - rayEx) / rayVx;
                        ip[1] = rayEy + t * rayVy;
                        ip[2] = rayEz + t * rayVz;

                        if (ip[1] < object.max[1] && ip[1] > object.min[1] && ip[2] < object.max[2] && ip[2] > object.min[2]) {
                                return true;
                        }
                }
                // top and bottom
                if (Math.abs(rayVy) > 1E-5) {
                        //top
                        t = (object.max[2] - rayEy) / rayVy;
                        ip[0] = rayEx + t * rayVx;
                        ip[2] = rayEz + t * rayVz;

                        if (ip[0] < object.max[0] && ip[0] > object.min[0] && ip[2] < object.max[2] && ip[2] > object.min[2]) {
                                return true;
                        }

                        //bottom
                        t = (object.min[2] - rayEy) / rayVy;
                        ip[0] = rayEx + t * rayVx;
                        ip[2] = rayEz + t * rayVz;

                        if (ip[0] < object.max[0] && ip[0] > object.min[0] && ip[2] < object.max[2] && ip[2] > object.min[2]) {
                                return true;
                        }
                }
                return false;
        }

        private void prepareMeshData() {
                RTFile scene;

                System.out.println("Vorverarbeitung 1 läuft");

                float[] p1, p2, p3;

                for (int objectsNumber = 0; objectsNumber < objects.size(); objectsNumber++) {
                        scene = objects.get(objectsNumber);

                        if (scene.getHeader().equals("TRIANGLE_MESH")) {
                                T_Mesh mesh = (T_Mesh) scene;

                                // init memory 
                                mesh.triangleNormals = new float[mesh.triangles.length][3];
                                mesh.triangleAreas = new float[mesh.triangles.length];

                                for (int i = 0; i < mesh.triangles.length; i++) {
                                        p1 = mesh.vertices[mesh.triangles[i][0]];
                                        p2 = mesh.vertices[mesh.triangles[i][1]];
                                        p3 = mesh.vertices[mesh.triangles[i][2]];
                                        // calculate and store triangle normal n and triangle area a
                                        mesh.triangleAreas[i] = calculateN(mesh.triangleNormals[i], p1, p2, p3);
                                }
                        }
                }
                System.out.println("Vorverarbeitung 1 beendet");
        }

        // view dependend precalculation dependend on type of mesh shading
        // vertexNormals for phong and gouraud shading
        // vertexColors for gouraud shading
        // triangleColors for flat lighting
        private void precalculateMeshDataShading() {
                RTFile scene;

                System.out.println("Vorverarbeitung 2 läuft");

                float rayEx, rayEy, rayEz, rayVx, rayVy, rayVz;
                double rayVn;
                Color color;
                float x, y, z;
                float[] ip = new float[3];
                float[] n = new float[3];
                float[] l = new float[3];
                float[] v = new float[3];
                float[] material;
                float materialN;
                int matIndex;

                for (int objectsNumber = 0; objectsNumber < objects.size(); objectsNumber++) {
                        scene = objects.get(objectsNumber);

                        if (scene.getHeader().equals("TRIANGLE_MESH")) {
                                T_Mesh mesh = (T_Mesh) scene;

                                switch (mesh.fgp) {
                                        case 'f':
                                        case 'F':
                                                // for flat-shading: initialize and calculate triangle
                                                // colors
                                                mesh.triangleColors = new float[mesh.triangles.length][3];

                                                rayEx = 0.0f;
                                                rayEy = 0.0f;
                                                rayEz = 0.0f;

                                                // loop over all triangles
                                                for (int i = 0; i < mesh.triangles.length; i++) {
                                                        // the intersection point is the first vertex of the
                                                        // triangle
                                                        ip = mesh.vertices[mesh.triangles[i][0]];

                                                        // the material is the material of the first triangle
                                                        // point
                                                        matIndex = mesh.verticesMat[mesh.triangles[i][0]];
                                                        material = mesh.materials[matIndex];
                                                        materialN = mesh.materialsN[matIndex];

                                                        // x, y, z: view coordinates are intersection point
                                                        x = ip[0];
                                                        y = ip[1];
                                                        z = ip[2];

                                                        // ray vector
                                                        rayVx = x - rayEx;
                                                        rayVy = y - rayEy;
                                                        rayVz = z - rayEz;

                                                        // fetch precalculated face normal
                                                        n = mesh.triangleNormals[i];

                                                        rayVn = rayVx * n[0] + rayVy * n[1] + rayVz * n[2];

                                                        // backface? => next triangle
                                                        if (rayVn >= 0) {
                                                                continue;
                                                        }

                                                        // light vector at the intersection point
                                                        l[0] = ICenter[0] - ip[0];
                                                        l[1] = ICenter[1] - ip[1];
                                                        l[2] = ICenter[2] - ip[2];
                                                        normalize(l);

                                                        // viewing vector at intersection point
                                                        v[0] = -rayVx;
                                                        v[1] = -rayVy;
                                                        v[2] = -rayVz;
                                                        normalize(v);

                                                        // illuminate
                                                        color = phongIlluminate(material, materialN, l, n, v, Ia, Ids);

                                                        // write color to triangle
                                                        mesh.triangleColors[i][0] = color.getRed() / 255.0f;
                                                        mesh.triangleColors[i][1] = color.getGreen() / 255.0f;
                                                        mesh.triangleColors[i][2] = color.getBlue() / 255.0f;
                                                }

                                                break;

                                        case 'p':
                                        case 'P':
                                        case 'g':
                                        case 'G':
                                                // initialize and calculate averaged vertex normals
                                                mesh.vertexNormals = new float[mesh.vertices.length][3];

                                                // loop over all vertices to initialize
                                                for (int j = 0; j < mesh.vertices.length; j++) {
                                                        for (int k = 0; k < 3; k++) {
                                                                mesh.vertexNormals[j][k] = 0.0f;
                                                        }
                                                }

                                                // loop over all faces to contribute
                                                for (int i = 0; i < mesh.triangles.length; i++) {
                                                        for (int j = 0; j < 3; j++) {
                                                                for (int k = 0; k < 3; k++) // loop over all vertices to normalize
                                                                {
                                                                        mesh.vertexNormals[mesh.triangles[i][j]][k] += mesh.triangleNormals[i][k];
                                                                }
                                                        }
                                                }
                                                
                                                for (int j = 0; j < mesh.vertices.length; j++) {
                                                        normalize(mesh.vertexNormals[j]);
                                                }

                                                // these are all preparations for phong shading
                                                if (mesh.fgp == 'p' || mesh.fgp == 'P') {
                                                        break;
                                                }

                                                // for gouraud-shading: initialize and calculate vertex
                                                // colors
                                                mesh.vertexColors = new float[mesh.vertices.length][3];

                                                rayEx = 0.0f;
                                                rayEy = 0.0f;
                                                rayEz = 0.0f;

                                                // loop over all vertices
                                                for (int i = 0; i < mesh.vertices.length; i++) {
                                                        // the intersection point is the vertex
                                                        ip = mesh.vertices[i];

                                                        // the material is the material of the vertex
                                                        matIndex = mesh.verticesMat[i];
                                                        material = mesh.materials[matIndex];
                                                        materialN = mesh.materialsN[matIndex];

                                                        // x, y, z: view coordinates are intersection point
                                                        x = ip[0];
                                                        y = ip[1];
                                                        z = ip[2];

                                                        // ray vector
                                                        rayVx = x - rayEx;
                                                        rayVy = y - rayEy;
                                                        rayVz = z - rayEz;

                                                        // fetch precalculated vertex normal
                                                        n = mesh.vertexNormals[i];

                                                        rayVn = rayVx * n[0] + rayVy * n[1] + rayVz * n[2];

                                                        // backface? => next vertex
                                                        if (rayVn >= 0) {
                                                                continue;
                                                        }

                                                        // light vector at the intersection point
                                                        l[0] = ICenter[0] - ip[0];
                                                        l[1] = ICenter[1] - ip[1];
                                                        l[2] = ICenter[2] - ip[2];
                                                        normalize(l);

                                                        // viewing vector at intersection point
                                                        v[0] = -rayVx;
                                                        v[1] = -rayVy;
                                                        v[2] = -rayVz;
                                                        normalize(v);

                                                        // illuminate
                                                        color = phongIlluminate(material, materialN, l, n, v, Ia, Ids);

                                                        // write color to vertex
                                                        mesh.vertexColors[i][0] = color.getRed() / 255.0f;
                                                        mesh.vertexColors[i][1] = color.getGreen() / 255.0f;
                                                        mesh.vertexColors[i][2] = color.getBlue() / 255.0f;
                                                }
                                }
                        }
                }
                System.out.println("Vorverarbeitung 2 beendet");
        }

        private float calculateN(float[] fn, float[] p1, float[] p2, float[] p3) {

                float vecAx, vecAy, vecAz;
                vecAx = p2[0] - p1[0];
                vecAy = p2[1] - p1[1];
                vecAz = p2[2] - p1[2];

                float vecBx, vecBy, vecBz;
                vecBx = p3[0] - p1[0];
                vecBy = p3[1] - p1[1];
                vecBz = p3[2] - p1[2];

                //crossproduct
                fn[0] = vecAy * vecBz - vecAz * vecBy;
                fn[1] = -(vecAx * vecBz - vecAz * vecBx);
                fn[2] = vecAx * vecBy - vecAy * vecBx;

                return normalize(fn) / 2;
        }

        private void calculateN2(float[] fn, float[] p1, float[] p2, float[] p3) {

                float vecAx, vecAy, vecAz;
                vecAx = p2[0] - p1[0];
                vecAy = p2[1] - p1[1];
                vecAz = p2[2] - p1[2];

                float vecBx, vecBy, vecBz;
                vecBx = p3[0] - p1[0];
                vecBy = p3[1] - p1[1];
                vecBz = p3[2] - p1[2];

                //crossproduct
                fn[0] = vecAy * vecBz - vecAz * vecBy;
                fn[1] = -(vecAx * vecBz - vecAz * vecBx);
                fn[2] = vecAx * vecBy - vecAy * vecBx;

        }

        public static void main(String[] args) {
                Raytracer00 rt = new Raytracer00();

                //rt.doRayTrace();
        }
}
