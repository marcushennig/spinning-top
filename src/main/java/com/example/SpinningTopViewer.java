package com.example;

import com.jme3.app.SimpleApplication;
import com.jme3.light.AmbientLight;
import com.jme3.light.DirectionalLight;
import com.jme3.material.Material;
import com.jme3.math.ColorRGBA;
import com.jme3.font.BitmapText;
import com.jme3.math.Quaternion;
import com.jme3.math.Vector3f;
import com.jme3.scene.Geometry;
import com.jme3.scene.debug.Arrow;
import com.jme3.scene.Mesh;
import com.jme3.scene.VertexBuffer;
import com.jme3.util.BufferUtils;
import com.jme3.system.AppSettings;

/**
 * Visualises the spinning top using jMonkeyEngine.
 */
public class SpinningTopViewer extends SimpleApplication {

    private final SpinningTopSimulation simulation = new SpinningTopSimulation();
    private double time = 0.0;
    private Geometry topGeom;
    private BitmapText hudText;
    private Geometry axisArrowGeom;
    private Geometry trajectoryGeom;
    private final java.util.List<Vector3f> trajectory = new java.util.ArrayList<>();

    public static void main(String[] args) {
        var app = new SpinningTopViewer();
        var settings = new AppSettings(true);
        
        settings.setTitle("Spinning Top Simulation"); // Set window title here
        settings.setSamples(16); // 4x MSAA (anti-aliasing)
        
        app.setSettings(settings);
        app.start();
    }

    @Override
    public void simpleInitApp() {
        var surface = new DupineCyclideSurface(30f, 27.8f, 5.99f, 9f, 100, 100);
        var mesh = surface.createMesh();
        topGeom = new Geometry("top", mesh);

        viewPort.setBackgroundColor(ColorRGBA.White);
        setDisplayStatView(false); // Hides FPS and stats overlay

        var mat = new Material(assetManager, "Common/MatDefs/Light/Lighting.j3md");
        mat.setBoolean("UseMaterialColors", true);
        mat.setColor("Diffuse", ColorRGBA.Orange);     // or interpolate via vertex color
        mat.setColor("Specular", ColorRGBA.White);     // white specular highlights
        mat.setFloat("Shininess", 64f);

        topGeom.setMaterial(mat);
        rootNode.attachChild(topGeom);
 
        // Directional Light (Sun-like light)
        var sun = new DirectionalLight();
        sun.setColor(ColorRGBA.White);
        sun.setDirection(new Vector3f(-1, -2, -3).normalizeLocal());  // light from above-left
        rootNode.addLight(sun);
        // Ambient Light (Soft fill)
        var ambient = new AmbientLight();
        ambient.setColor(ColorRGBA.White.mult(0.2f));  // subtle base illumination
        rootNode.addLight(ambient);

        
        cam.setLocation(new Vector3f(0, 0, 100));
        cam.lookAt(Vector3f.ZERO, Vector3f.UNIT_Y);

        var font = assetManager.loadFont("Interface/Fonts/Default.fnt");
        hudText = new BitmapText(font);
        hudText.setSize(font.getCharSet().getRenderedSize());
        hudText.setColor(ColorRGBA.Black);
        hudText.setLocalTranslation(10f, cam.getHeight() - 10f, 0);
        guiNode.attachChild(hudText);

        // Red arrow indicating the current (x,y,z) position of the centre of mass
        var arrow = new Arrow(Vector3f.UNIT_Z.mult((float) simulation.getConstants().l));
        arrow.setLineWidth(2);
        axisArrowGeom = new Geometry("comArrow", arrow);
        var redMat = new Material(assetManager, "Common/MatDefs/Misc/Unshaded.j3md");
        redMat.setColor("Color", ColorRGBA.Red);
        axisArrowGeom.setMaterial(redMat);
        rootNode.attachChild(axisArrowGeom);

        // Blue polyline for the recent trajectory
        trajectoryGeom = new Geometry("trajectory", new Mesh());
        var blueMat = new Material(assetManager, "Common/MatDefs/Misc/Unshaded.j3md");
        blueMat.setColor("Color", ColorRGBA.Blue);
        trajectoryGeom.setMaterial(blueMat);
        rootNode.attachChild(trajectoryGeom);
    }
    /**
     * Resets the simulation to initial conditions.
     * This method is called when the application starts.
     * @param tpf time per frame, not used in this method
     */
    @Override
    public void simpleUpdate(float tpf) {
        
        double dt = simulation.getDt(); // Integration time step
        simulation.evolute(time, time + dt);
        time += dt;

        simulation.updateEnergy();

        var s = simulation.getState();
        var rot = new Quaternion().fromAngleAxis((float) s.phi, Vector3f.UNIT_Y)
                .mult(new Quaternion().fromAngleAxis((float) s.theta, Vector3f.UNIT_X))
                .mult(new Quaternion().fromAngleAxis((float) s.psi, Vector3f.UNIT_Z));
        topGeom.setLocalRotation(rot);

        // Compute centre-of-mass position from Euler angles
        var comPos = rot.mult(new Vector3f(0f, 0f, (float) simulation.getConstants().l));
        ((Arrow) axisArrowGeom.getMesh()).setArrowExtent(comPos);
        axisArrowGeom.updateModelBound();

        updateTrajectory(comPos);
        updateTrajectoryGeometry();

        hudText.setText(String.format(
                "Time: %.3f\nEnergy: %.3f\ntheta: %.3f\nphi: %.3f\npsi: %.3f",
                time, s.energy, s.theta, s.phi, s.psi));
    }

    private void updateTrajectory(Vector3f pos) {
        trajectory.add(pos.clone());
        if (trajectory.size() > 100) {
            trajectory.remove(0);
        }
    }

    private void updateTrajectoryGeometry() {
        var mesh = new Mesh();
        mesh.setMode(Mesh.Mode.LineStrip);
        var pts = trajectory.toArray(Vector3f[]::new);
        if (pts.length < 2) {
            trajectoryGeom.setMesh(mesh);
            return;
        }
        int[] indices = new int[pts.length];
        for (int i = 0; i < pts.length; i++) {
            indices[i] = i;
        }
        mesh.setBuffer(VertexBuffer.Type.Position, 3, BufferUtils.createFloatBuffer(pts));
        mesh.setBuffer(VertexBuffer.Type.Index, 1, BufferUtils.createIntBuffer(indices));
        mesh.updateBound();
        mesh.updateCounts();
        trajectoryGeom.setMesh(mesh);
    }

}
