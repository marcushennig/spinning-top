package com.example;

import org.lwjgl.*;
import org.lwjgl.glfw.*;
import org.lwjgl.opengl.*;
import org.lwjgl.system.*;

import java.nio.*;

import static org.lwjgl.glfw.Callbacks.*;
import static org.lwjgl.glfw.GLFW.*;
import static org.lwjgl.opengl.GL11.*;
import static org.lwjgl.system.MemoryStack.*;
import static org.lwjgl.system.MemoryUtil.*;

import org.lwjgl.opengl.GL;
import java.util.logging.Logger;

public class SpinningTopViewer {

    private static final Logger logger = Logger.getLogger(SpinningTopViewer.class.getName());

    private SpinningTopSimulation simulation = new SpinningTopSimulation();
    private double time = 0;
    private long window;
    private final int width = 400;
    private final int height = 400;

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
        //
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
            glRotated(simulation.getParameter().phi * 180 / Math.PI, 0, 1, 0);
            glRotated(simulation.getParameter().theta * 180 / Math.PI, 1, 0, 0);
            glRotated(simulation.getParameter().psi * 180 / Math.PI, 0, 0, 1);

            renderTorus(20, 10, 100, 100);

            time += simulation.getParameter().dt;
            simulation.NDSolve(time, time + simulation.getParameter().dt);

            glfwSwapBuffers(window); // swap the color buffers
            // Poll for window events. The key callback above will only be
			// invoked during this call.
			glfwPollEvents();
        }
        logger.info("Exited main loop.");
    }

    /**
     * Renders a torus with the specified major and minor radii,
     * and the number of segments for each.
     * @param majorRadius The radius of the torus' center circle.
     * @param minorRadius The radius of the torus' tube.
     * @param numMajor The number of segments around the major circle.
     * @param numMinor The number of segments around the minor circle.
     */
    private void renderTorus(
        float majorRadius, 
        float minorRadius, 
        int numMajor, 
        int numMinor){
        
        float majorStep = (float)(2.0f * Math.PI / numMajor);
        float minorStep = (float)(2.0f * Math.PI / numMinor);

        glPushMatrix();

        // Enable solid fill
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        // Set material color
        glColor3f(1.0f, 0.6f, 0.2f);

        // Loop through torus surface
        for (int i = 0; i < numMajor; ++i) {
            float a0 = i * majorStep;
            float a1 = (i + 1) * majorStep;
            float cosA0 = (float) Math.cos(a0);
            float sinA0 = (float) Math.sin(a0);
            float cosA1 = (float) Math.cos(a1);
            float sinA1 = (float) Math.sin(a1);

            glBegin(GL_TRIANGLE_STRIP);
            for (int j = 0; j <= numMinor; ++j) {
                float b = j * minorStep;
                float cosB = (float) Math.cos(b);
                float sinB = (float) Math.sin(b);

                float r = minorRadius;

                // First vertex
                float nz = cosA0 * cosB;
                float nx = sinA0 * cosB;
                float ny = sinB;
                float z = (majorRadius + r * cosB) * cosA0;
                float x = (majorRadius + r * cosB) * sinA0;
                float y = r * sinB;
                glNormal3f(nx, ny, nz);
                glVertex3f(x, y, z);

                // Second vertex
                nz = cosA1 * cosB;
                nx = sinA1 * cosB;
                ny = sinB;
                z = (majorRadius + r * cosB) * cosA1;
                x = (majorRadius + r * cosB) * sinA1;
                y = r * sinB;
                glNormal3f(nx, ny, nz);
                glVertex3f(x, y, z);
            }
            glEnd();
        }
        glPopMatrix();
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
