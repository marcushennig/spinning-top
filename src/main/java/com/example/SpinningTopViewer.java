package com.example;

import com.jme3.app.SimpleApplication;
import com.jme3.light.DirectionalLight;
import com.jme3.material.Material;
import com.jme3.math.ColorRGBA;
import com.jme3.math.FastMath;
import com.jme3.math.Quaternion;
import com.jme3.math.Vector3f;
import com.jme3.scene.Geometry;
import com.jme3.scene.Mesh;
import com.jme3.scene.VertexBuffer;
import com.jme3.util.BufferUtils;

/**
 * Visualises the spinning top using jMonkeyEngine.
 */
public class SpinningTopViewer extends SimpleApplication {

    private final SpinningTopSimulation simulation = new SpinningTopSimulation();
    private double time = 0.0;
    private Geometry topGeom;

    public static void main(String[] args) {
        SpinningTopViewer app = new SpinningTopViewer();
        app.start();
    }

    @Override
    public void simpleInitApp() {
        Mesh mesh = createDupinCyclide(30f, 27.8f, 5.99f, 9f, 100, 100);
        topGeom = new Geometry("top", mesh);

        Material mat = new Material(assetManager, "Common/MatDefs/Light/Lighting.j3md");
        mat.setBoolean("UseMaterialColors", true);
        mat.setColor("Diffuse", new ColorRGBA(1f, 0.6f, 0.2f, 1f));
        mat.setColor("Ambient", ColorRGBA.White);
        topGeom.setMaterial(mat);
        rootNode.attachChild(topGeom);

        DirectionalLight light = new DirectionalLight();
        light.setDirection(new Vector3f(-1, -1, -1).normalizeLocal());
        rootNode.addLight(light);

        cam.setLocation(new Vector3f(0, 0, 100));
        cam.lookAt(Vector3f.ZERO, Vector3f.UNIT_Y);
    }

    @Override
    public void simpleUpdate(float tpf) {
        simulation.evolute(time, time + tpf);
        time += tpf;

        SpinningTopState s = simulation.getState();
        Quaternion rot = new Quaternion().fromAngleAxis((float) s.phi, Vector3f.UNIT_Y)
                .mult(new Quaternion().fromAngleAxis((float) s.theta, Vector3f.UNIT_X))
                .mult(new Quaternion().fromAngleAxis((float) s.psi, Vector3f.UNIT_Z));
        topGeom.setLocalRotation(rot);
    }

    /**
     * Build a Dupin cyclide mesh similar to the old OpenGL renderer.
     */
    private Mesh createDupinCyclide(float a, float b, float c, float d, int numU, int numV) {
        int vertCount = (numU + 1) * (numV + 1);
        Vector3f[] vertices = new Vector3f[vertCount];
        Vector3f[] normals = new Vector3f[vertCount];

        for (int i = 0; i <= numU; i++) {
            float u = i * FastMath.TWO_PI / numU;
            for (int j = 0; j <= numV; j++) {
                float v = j * FastMath.TWO_PI / numV;
                int index = i * (numV + 1) + j;
                vertices[index] = dupinCyclidePoint(a, b, c, d, u, v);
                normals[index] = dupinCyclideNormal(a, b, c, d, u, v);
            }
        }

        int[] indices = new int[numU * numV * 6];
        int idx = 0;
        for (int i = 0; i < numU; i++) {
            for (int j = 0; j < numV; j++) {
                int p0 = i * (numV + 1) + j;
                int p1 = (i + 1) * (numV + 1) + j;
                int p2 = (i + 1) * (numV + 1) + j + 1;
                int p3 = i * (numV + 1) + j + 1;
                indices[idx++] = p0;
                indices[idx++] = p1;
                indices[idx++] = p2;
                indices[idx++] = p0;
                indices[idx++] = p2;
                indices[idx++] = p3;
            }
        }

        Mesh mesh = new Mesh();
        mesh.setBuffer(VertexBuffer.Type.Position, 3, BufferUtils.createFloatBuffer(vertices));
        mesh.setBuffer(VertexBuffer.Type.Normal, 3, BufferUtils.createFloatBuffer(normals));
        mesh.setBuffer(VertexBuffer.Type.Index, 3, BufferUtils.createIntBuffer(indices));
        mesh.updateBound();
        return mesh;
    }

    private Vector3f dupinCyclidePoint(float a, float b, float c, float d, float u, float v) {
        float cu = FastMath.cos(u);
        float su = FastMath.sin(u);
        float cv = FastMath.cos(v);
        float sv = FastMath.sin(v);
        float denom = a - c * cu * cv;
        float x = (d * (c - a * cu * cv) + b * b * cu) / denom;
        float y = b * su * (a - d * cv) / denom;
        float z = b * sv * (c * cu - d) / denom;
        return new Vector3f(x, y, z);
    }

    private Vector3f dupinCyclideNormal(float a, float b, float c, float d, float u, float v) {
        float eps = 0.0001f;
        Vector3f p = dupinCyclidePoint(a, b, c, d, u, v);
        Vector3f pu = dupinCyclidePoint(a, b, c, d, u + eps, v);
        Vector3f pv = dupinCyclidePoint(a, b, c, d, u, v + eps);
        Vector3f du = pu.subtract(p);
        Vector3f dv = pv.subtract(p);
        return du.cross(dv).normalizeLocal();
    }
}
