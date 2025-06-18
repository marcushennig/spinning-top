package com.example;

/** Simple LWJGL viewer for the spinning top simulation. */
import static org.lwjgl.glfw.GLFW.*;
import static org.lwjgl.opengl.GL11.*;

import org.lwjgl.glfw.GLFWErrorCallback;
import org.lwjgl.opengl.GL;
import org.lwjgl.util.glu.Cylinder;
import org.lwjgl.util.glu.GLU;

public class SpinningTopViewer {
    private SpinningTopSimulation sim = new SpinningTopSimulation();
    private double time = 0;
    private long window;
    private final int width = 800;
    private final int height = 600;
    private final Cylinder cone = new Cylinder();
    private final GLU glu = new GLU();

    public static void main(String[] args) {
        new SpinningTopViewer().run();
    }

    public void run() {
        init();
        loop();
        glfwDestroyWindow(window);
        glfwTerminate();
    }

    private void init() {
        GLFWErrorCallback.createPrint(System.err).set();
        if (!glfwInit()) {
            throw new IllegalStateException("Unable to initialize GLFW");
        }
        window = glfwCreateWindow(width, height, "Spinning Top", 0, 0);
        if (window == 0) {
            throw new RuntimeException("Failed to create window");
        }
        glfwMakeContextCurrent(window);
        glfwSwapInterval(1);
        GL.createCapabilities();

        glClearColor(0f, 0f, 0f, 1f);
        glEnable(GL_DEPTH_TEST);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glu.gluPerspective(45.0f, (float) width / height, 1.0f, 500.0f);
        glMatrixMode(GL_MODELVIEW);
    }

    private void loop() {
        while (!glfwWindowShouldClose(window)) {
            glfwPollEvents();
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            glLoadIdentity();
            glu.gluLookAt(0, 0, 100, 0, 0, 0, 0, 1, 0);

            glRotated(sim.getParameter().phi * 180 / Math.PI, 0, 1, 0);
            glRotated(sim.getParameter().theta * 180 / Math.PI, 1, 0, 0);
            glRotated(sim.getParameter().psi * 180 / Math.PI, 0, 0, 1);
            drawTop();

            time += sim.getParameter().dt;
            sim.NDSolve(time, time + sim.getParameter().dt);

            glfwSwapBuffers(window);
        }
    }

    private void drawTop() {
        glColor3f(0.8f, 0.2f, 0.2f);
        cone.draw(5f, 0f, 20f, 20, 20);
    }
}
