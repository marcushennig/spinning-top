package com.example;

import com.jme3.app.SimpleApplication;
import com.jme3.light.DirectionalLight;
import com.jme3.material.Material;
import com.jme3.math.ColorRGBA;
import com.jme3.math.Quaternion;
import com.jme3.math.Vector3f;
import com.jme3.scene.Geometry;
import com.jme3.system.AppSettings;

/**
 * Visualises the spinning top using jMonkeyEngine.
 */
public class SpinningTopViewer extends SimpleApplication {

    private final SpinningTopSimulation simulation = new SpinningTopSimulation();
    private double time = 0.0;
    private Geometry topGeom;

    public static void main(String[] args) {
        var app = new SpinningTopViewer();
        var settings = new AppSettings(true);
        settings.setSamples(4); // 4x MSAA (anti-aliasing)
        app.setSettings(settings);
        app.start();
    }

    @Override
    public void simpleInitApp() {
        var surface = new DupineCyclideSurface(30f, 27.8f, 5.99f, 9f, 100, 100);
        var mesh = surface.createMesh();
        topGeom = new Geometry("top", mesh);

        viewPort.setBackgroundColor(ColorRGBA.White);

        var mat = new Material(assetManager, "Common/MatDefs/Light/Lighting.j3md");
        mat.setBoolean("UseMaterialColors", true);
        mat.setColor("Diffuse", new ColorRGBA(1f, 0.6f, 0.2f, 1f));
        mat.setColor("Ambient", ColorRGBA.White);
        // Make it shiny:
        mat.setColor("Specular", ColorRGBA.White); // Strong specular highlight
        mat.setFloat("Shininess", 64f); // Higher value = smaller, sharper highlight (default is 1–128)

        topGeom.setMaterial(mat);
        rootNode.attachChild(topGeom);

        var light = new DirectionalLight();
        light.setDirection(new Vector3f(-1, -1, -1).normalizeLocal());
        rootNode.addLight(light);

        cam.setLocation(new Vector3f(0, 0, 100));
        cam.lookAt(Vector3f.ZERO, Vector3f.UNIT_Y);
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

        var s = simulation.getState();
        var rot = new Quaternion().fromAngleAxis((float) s.phi, Vector3f.UNIT_Y)
                .mult(new Quaternion().fromAngleAxis((float) s.theta, Vector3f.UNIT_X))
                .mult(new Quaternion().fromAngleAxis((float) s.psi, Vector3f.UNIT_Z));
        topGeom.setLocalRotation(rot);
    }

}
