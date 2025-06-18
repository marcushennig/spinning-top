package com.example;

/** Simple LWJGL viewer for the spinning top simulation. */
import static org.lwjgl.glfw.GLFW.*;
import static org.lwjgl.opengl.GL11.*;

import org.lwjgl.glfw.GLFWErrorCallback;
import org.lwjgl.opengl.GL;

public class SpinningTopViewer {
    
    private SpinningTopSimulation sim = new SpinningTopSimulation();
    private double time = 0;
    private long window;
    private final int width = 800;
    private final int height = 600;

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
        perspective(45.0f, (float) width / height, 1.0f, 500.0f);
        glMatrixMode(GL_MODELVIEW);
    }

    private void loop() {
        while (!glfwWindowShouldClose(window)) {
            glfwPollEvents();
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            glLoadIdentity();
            lookAt(0, 0, 100, 0, 0, 0, 0, 1, 0);

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
    }

    private void perspective(float fovy, float aspect, float zNear, float zFar) {
        float fH = (float) Math.tan(Math.toRadians(fovy / 2)) * zNear;
        float fW = fH * aspect;
        glFrustum(-fW, fW, -fH, fH, zNear, zFar);
    }

    private void lookAt(
            double eyeX, double eyeY, double eyeZ,
            double centerX, double centerY, double centerZ,
            double upX, double upY, double upZ) {
        double[] forward = {
            centerX - eyeX,
            centerY - eyeY,
            centerZ - eyeZ
        };
        double fLen = Math.sqrt(forward[0]*forward[0] + forward[1]*forward[1] + forward[2]*forward[2]);
        forward[0] /= fLen; forward[1] /= fLen; forward[2] /= fLen;

        double[] up = { upX, upY, upZ };
        double upLen = Math.sqrt(up[0]*up[0] + up[1]*up[1] + up[2]*up[2]);
        up[0] /= upLen; up[1] /= upLen; up[2] /= upLen;

        double[] side = {
            forward[1]*up[2] - forward[2]*up[1],
            forward[2]*up[0] - forward[0]*up[2],
            forward[0]*up[1] - forward[1]*up[0]
        };
        double sLen = Math.sqrt(side[0]*side[0] + side[1]*side[1] + side[2]*side[2]);
        side[0] /= sLen; side[1] /= sLen; side[2] /= sLen;

        double[] up2 = {
            side[1]*forward[2] - side[2]*forward[1],
            side[2]*forward[0] - side[0]*forward[2],
            side[0]*forward[1] - side[1]*forward[0]
        };

        float[] m = new float[] {
            (float)side[0],   (float)up2[0],   (float)-forward[0],   0,
            (float)side[1],   (float)up2[1],   (float)-forward[1],   0,
            (float)side[2],   (float)up2[2],   (float)-forward[2],   0,
            0,                0,               0,                    1
        };
        glMultMatrixf(m);
        glTranslated(-eyeX, -eyeY, -eyeZ);
    }
}
