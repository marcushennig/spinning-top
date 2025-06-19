package com.example;

import com.jme3.math.FastMath;
import com.jme3.math.Vector3f;
import com.jme3.scene.Mesh;
import com.jme3.scene.VertexBuffer;
import com.jme3.util.BufferUtils;

/**
 * Represents a Dupin cyclide surface, which is a type of ruled surface
 * defined by four parameters (a, b, c, d) and discretized into a mesh.
 * This class provides methods to create the mesh and calculate points
 * and normals on the surface.
 * https://en.wikipedia.org/wiki/Dupin_cyclide
 */
public class DupineCyclideSurface {

    private final float a;
    private final float b;
    private final float c;
    private final float d;
    private final int numU;
    private final int numV;
    
    /**
     * Constructs a Dupin cyclide surface with the specified parameters.
     *
     * @param a    Parameter a of the cyclide.
     * @param b    Parameter b of the cyclide.
     * @param c    Parameter c of the cyclide.
     * @param d    Parameter d of the cyclide.
     * @param numU Number of segments in the u direction.
     * @param numV Number of segments in the v direction.
     */
    public DupineCyclideSurface(float a, float b, float c, float d, int numU, int numV) {
        this.a = a;
        this.b = b;
        this.c = c;
        this.d = d;
        this.numU = numU;
        this.numV = numV;
    }

    /**
     * Creates the mesh representing the Dupin cyclide.
     */
    public Mesh createMesh() {
        int vertCount = (numU + 1) * (numV + 1);
        Vector3f[] vertices = new Vector3f[vertCount];
        Vector3f[] normals = new Vector3f[vertCount];

        for (int i = 0; i <= numU; i++) {
            float u = i * FastMath.TWO_PI / numU;
            for (int j = 0; j <= numV; j++) {
                float v = j * FastMath.TWO_PI / numV;
                int index = i * (numV + 1) + j;
                vertices[index] = dupinCyclidePoint(u, v);
                normals[index] = dupinCyclideNormal(u, v);
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

    private Vector3f dupinCyclidePoint(float u, float v) {
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

    private Vector3f dupinCyclideNormal(float u, float v) {
        float eps = 0.0001f;
        Vector3f p = dupinCyclidePoint(u, v);
        Vector3f pu = dupinCyclidePoint(u + eps, v);
        Vector3f pv = dupinCyclidePoint(u, v + eps);
        Vector3f du = pu.subtract(p);
        Vector3f dv = pv.subtract(p);
        return du.cross(dv).normalizeLocal();
    }
}
