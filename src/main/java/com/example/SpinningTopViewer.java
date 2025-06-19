package com.example;

import org.lwjgl.*;
import org.lwjgl.glfw.*;
import org.lwjgl.system.*;

import java.nio.*;

import static org.lwjgl.glfw.Callbacks.*;
import static org.lwjgl.glfw.GLFW.*;
import static org.lwjgl.opengl.GL11.*;
import static org.lwjgl.system.MemoryStack.*;

import org.lwjgl.opengl.GL;
import java.util.logging.Logger;

public class SpinningTopViewer {

    private static final Logger logger = Logger.getLogger(SpinningTopViewer.class.getName());

    private SpinningTopSimulation simulation = new SpinningTopSimulation();
    private double time = 0;
    private long window;
    private final int width = 1000;
    private final int height = 600;

    /** Main entry point for the SpinningTopViewer application. */  
    public static void main(String[] args) {
        logger.info("Starting SpinningTopViewer...");
        new SpinningTopViewer().run();
        logger.info("SpinningTopViewer terminated.");
    }
    
    /** Main loop for the SpinningTopViewer. */
    public void run() {
        
        logger.info("Initializing...");
        init();
        
        logger.info("Entering main loop...");
        loop();
        
        logger.info("Destroying window and terminating GLFW...");
        glfwFreeCallbacks(window);
		glfwDestroyWindow(window);

		// Terminate GLFW and free the error callback
		glfwTerminate();
		glfwSetErrorCallback(null).free();
        
        logger.info("Exiting SpinningTopViewer.");
    }

    /** Initializes the OpenGL context and GLFW window. */
    private void init() {
        
        // Setup an error callback. The default implementation
        // will print the error message in System.err.
		GLFWErrorCallback.createPrint(System.err).set();

        // Initialize GLFW. Most GLFW functions will not work before doing this.
        if (!glfwInit()) {
            logger.severe("Unable to initialize GLFW");
            throw new IllegalStateException("Unable to initialize GLFW");
        }

        // Configure GLFW
		glfwDefaultWindowHints(); // optional, the current window hints are already the default
		glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE); // the window will stay hidden after creation
		glfwWindowHint(GLFW_RESIZABLE, GLFW_TRUE); // the window will be resizable

        
        // Create the window
        window = glfwCreateWindow(width, height, "Spinning Top", 0, 0);
        if (window == 0) {
            logger.severe("Failed to create window");
            throw new RuntimeException("Failed to create window");
        }
        
        // Setup a key callback. It will be called every time a key is pressed, repeated or released.
		glfwSetKeyCallback(window, (window, key, scancode, action, mods) -> {
			if ( key == GLFW_KEY_ESCAPE && action == GLFW_RELEASE )
				glfwSetWindowShouldClose(window, true); // We will detect this in the rendering loop
		});

        // Get the thread stack and push a new frame
		try ( MemoryStack stack = stackPush() ) {
			IntBuffer pWidth = stack.mallocInt(1); // int*
			IntBuffer pHeight = stack.mallocInt(1); // int*

			// Get the window size passed to glfwCreateWindow
			glfwGetWindowSize(window, pWidth, pHeight);

			// Get the resolution of the primary monitor
			GLFWVidMode vidmode = glfwGetVideoMode(glfwGetPrimaryMonitor());

			// Center the window
			glfwSetWindowPos(
				window,
				(vidmode.width() - pWidth.get(0)) / 2,
				(vidmode.height() - pHeight.get(0)) / 2
			);
		} // the stack frame is popped automatically

		// Make the OpenGL context current
		glfwMakeContextCurrent(window);
		logger.info("OpenGL context made current.");

        // Enable v-sync
		glfwSwapInterval(1);
        logger.info("OpenGL context made current and v-sync enabled.");
		
        // Make the window visible
		glfwShowWindow(window);
        logger.info("Window created and made visible.");
    }

    /** Main loop for rendering and updating the simulation. */
    private void loop() {

        // This line is critical for LWJGL's interoperation with GLFW's
		// OpenGL context, or any context that is managed externally.
		// LWJGL detects the context that is current in the current thread,
		// creates the GLCapabilities instance and makes the OpenGL
		// bindings available for use.
		GL.createCapabilities();

        // Set the clear color
		glClearColor(1f, 1f, 1f, 1f);
   
        glEnable(GL_DEPTH_TEST);     // Enables Z-buffering
        glEnable(GL_LIGHTING);       // Enables lighting
        glEnable(GL_LIGHT0);         // Enables a default light
        glEnable(GL_NORMALIZE);      // Normalize normals after scaling
        glShadeModel(GL_SMOOTH);    // Smooth shading

        // Set up light properties
        // Ambient light
        // Diffuse light
        // Position of the light source
        FloatBuffer lightAmbient = BufferUtils.createFloatBuffer(4).put(new float[]{0.2f, 0.2f, 0.2f, 1.0f});
        FloatBuffer lightDiffuse = BufferUtils.createFloatBuffer(4).put(new float[]{0.8f, 0.8f, 0.8f, 1.0f});
        FloatBuffer lightPosition = BufferUtils.createFloatBuffer(4).put(new float[]{1.0f, 1.0f, 1.0f, 0.0f}); // directional light

        lightAmbient.flip();
        lightDiffuse.flip();
        lightPosition.flip();

        glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmbient);
        glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiffuse);
        glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);

        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
        glEnable(GL_COLOR_MATERIAL);
        FloatBuffer matSpecular = BufferUtils.createFloatBuffer(4).put(new float[]{1.0f, 1.0f, 1.0f, 1.0f});
        matSpecular.flip();
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, matSpecular);
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 64.0f); 

        glEnable(GL_DEPTH_TEST);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        perspective(45.0f, (float) width / height, 1.0f, 500.0f);
        glMatrixMode(GL_MODELVIEW);

        // Run the rendering loop until the user has attempted to close
		// the window or has pressed the ESCAPE key.
        while ( !glfwWindowShouldClose(window) ) {
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // clear the framebuffer

            glLoadIdentity();
            
            lookAt(0, 0, 100, 0, 0, 0, 0, 1, 0);
            glRotated(simulation.getState().phi * 180 / Math.PI, 0, 1, 0);
            glRotated(simulation.getState().theta * 180 / Math.PI, 1, 0, 0);
            glRotated(simulation.getState().psi * 180 / Math.PI, 0, 0, 1);

            renderDupinCyclide(30.0f, 27.8f, 5.99f, 9f, 100, 100);

            time += simulation.getDt();
            simulation.NDSolve(time, time + simulation.getDt());

            glfwSwapBuffers(window); // swap the color buffers
            // Poll for window events. The key callback above will only be
			// invoked during this call.
			glfwPollEvents();
        }
        logger.info("Exited main loop.");
    }
    /**
     * Renders a Dupin cyclide using triangle strips.
     * A Dupin cyclide is a type of surface generated by the rotation of a circle
     * around a fixed axis, where the circle's radius varies along the axis.
     * @param a the radius of the first circle
     * @param b the radius of the second circle
     * @param c the radius of the third circle 
     * @param d the radius of the fourth circle
     * @param numU the number of segments in the u direction
     * @param numV  the number of segments in the v direction
     */
    private void renderDupinCyclide(
        float a, float b, float c, float d,
        int numU, int numV) {

        float uStep = (float)(2.0f * Math.PI / numU);
        float vStep = (float)(2.0f * Math.PI / numV);

        glPushMatrix();
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        // Set material color
        glColor3f(1.0f, 0.6f, 0.2f);

        for (int i = 0; i < numU; ++i) {
            float u0 = i * uStep;
            float u1 = (i + 1) * uStep;

            glBegin(GL_TRIANGLE_STRIP);
            for (int j = 0; j <= numV; ++j) {
                float v = j * vStep;

                float[] p1 = dupinCyclidePoint(a, b, c, d, u0, v);
                float[] n1 = dupinCyclideNormal(a, b, c, d, u0, v);
                glNormal3f(n1[0], n1[1], n1[2]);
                glVertex3f(p1[0], p1[1], p1[2]);

                float[] p2 = dupinCyclidePoint(a, b, c, d, u1, v);
                float[] n2 = dupinCyclideNormal(a, b, c, d, u1, v);
                glNormal3f(n2[0], n2[1], n2[2]);
                glVertex3f(p2[0], p2[1], p2[2]);
            }
            glEnd();
        }

        glPopMatrix();
    }
    
    /**
     * Calculates the normal vector for a point on a Dupin cyclide surface.
     * The Dupin cyclide is defined by four parameters (a, b, c, d) and two parameters (u, v).
     * The normal vector is computed using the partial derivatives of the surface
     * with respect to u and v, and then normalizing the resulting
     * @param a the radius of the first circle
     * @param b the radius of the second circle
     * @param c the radius of the third circle
     * @param d the radius of the fourth circle
     * @param u the parameter in the u direction (0 to 2π)
     * @param v the parameter in the v direction (0 to 2π)
     * @return the normal vector at the point (u, v) on the Dupin cyclide surface
     */
    private float[] dupinCyclideNormal(float a, float b, float c, float d, float u, float v) {
        float eps = 0.0001f;

        float[] p = dupinCyclidePoint(a, b, c, d, u, v);
        float[] pu = dupinCyclidePoint(a, b, c, d, u + eps, v);
        float[] pv = dupinCyclidePoint(a, b, c, d, u, v + eps);

        float[] du = new float[]{pu[0] - p[0], pu[1] - p[1], pu[2] - p[2]};
        float[] dv = new float[]{pv[0] - p[0], pv[1] - p[1], pv[2] - p[2]};

        float[] normal = cross(du, dv);
        normalize(normal);
        return normal;
    }

    /**
     * Calculates the 3D coordinates of a point on a Dupin cyclide surface.
     * @param a the radius of the first circle
     * @param b the radius of the second circle
     * @param c the radius of the third circle
     * @param d the radius of the fourth circle
     * @param u the parameter in the u direction (0 to 2π)
     * @param v the parameter in the v direction (0 to 2π)
     * @return the 3D coordinates of the point on the Dupin cyclide surface
     */
    private float[] dupinCyclidePoint(float a, float b, float c, float d, float u, float v) {
        
        float cu = (float) Math.cos(u);
        float su = (float) Math.sin(u);
        float cv = (float) Math.cos(v);
        float sv = (float) Math.sin(v);

        float denom = a - c * cu * cv;

        float x = (d * (c - a * cu * cv) + b * b * cu) / denom;
        float y = b * su * (a - d * cv) / denom;
        float z = b * sv * (c * cu - d) / denom;

        return new float[]{x, y, z};
    }

    /**
     * Computes the cross product of two 3D vectors.
     * This method takes two 3D vectors represented as float arrays and returns their cross product
     * @param a the first vector
     * @param b the second vector
     * @return the cross product of a and b as a float array
     */
    private float[] cross(float[] a, float[] b) {
        return new float[] {
            a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0]
        };
    }

    /**
     * Normalizes a 3D vector represented as a float array.
     * This method modifies the input vector to have a length of 1.
     * @param v the vector to normalize
     */
    private void normalize(float[] v) {
        float length = (float) Math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        if (length != 0.0f) {
            v[0] /= length;
            v[1] /= length;
            v[2] /= length;
        }
    }

    /**
     *  Sets up a perspective projection matrix.
     *  This method calculates the field of view, aspect ratio, and near/far clipping planes,
     *  and applies the perspective projection to the OpenGL context.
     * @param fovy the field of view in the y direction, in degrees
     * @param aspect the aspect ratio of the viewport (width / height)
     * @param zNear the distance to the near clipping plane
     * @param zFar the distance to the far clipping plane
     */
    private void perspective(float fovy, float aspect, float zNear, float zFar) {
        float fH = (float) Math.tan(Math.toRadians(fovy / 2)) * zNear;
        float fW = fH * aspect;
        glFrustum(-fW, fW, -fH, fH, zNear, zFar);
    }
    
    /**
     * Sets the camera to look at a specific point in 3D space.
     * This method computes the forward, up, and side vectors based on the
     * eye position, center position, and up vector, and applies the view
     * @param eyeX   
     * @param eyeY
     * @param eyeZ
     * @param centerX
     * @param centerY
     * @param centerZ
     * @param upX
     * @param upY
     * @param upZ
     */
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
