package com.example;

/** Simple JOGL viewer for the spinning top simulation. */
import com.jogamp.opengl.*;
import com.jogamp.opengl.glu.GLU;
import com.jogamp.opengl.util.gl2.GLUT;
import com.jogamp.opengl.awt.GLCanvas;
import com.jogamp.opengl.util.Animator;
import javax.swing.JFrame;

public class SpinningTopViewer implements GLEventListener {
    private SpinningTopSimulation sim = new SpinningTopSimulation();
    private double time = 0;

    public static void main(String[] args){
        GLProfile profile = GLProfile.get(GLProfile.GL2);
        GLCapabilities caps = new GLCapabilities(profile);
        GLCanvas canvas = new GLCanvas(caps);
        SpinningTopViewer app = new SpinningTopViewer();
        canvas.addGLEventListener(app);
        Animator animator = new Animator(canvas);
        JFrame frame = new JFrame("Spinning Top");
        frame.getContentPane().add(canvas);
        frame.setSize(800,600);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setVisible(true);
        animator.start();
    }

    @Override
    public void init(GLAutoDrawable drawable){
        GL2 gl = drawable.getGL().getGL2();
        gl.glClearColor(0,0,0,1);
        gl.glEnable(GL.GL_DEPTH_TEST);
    }

    @Override
    public void reshape(GLAutoDrawable drawable,int x,int y,int width,int height){
        GL2 gl = drawable.getGL().getGL2();
        gl.glViewport(0,0,width,height);
        gl.glMatrixMode(GL2.GL_PROJECTION);
        gl.glLoadIdentity();
        float aspect = (float)width/(float)height;
        glu.gluPerspective(45.0f, aspect, 1.0f, 500.0f);
        gl.glMatrixMode(GL2.GL_MODELVIEW);
    }

    private GLU glu = new GLU();

    @Override
    public void display(GLAutoDrawable drawable){
        GL2 gl = drawable.getGL().getGL2();
        gl.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT);
        gl.glLoadIdentity();
        glu.gluLookAt(0,0,100, 0,0,0, 0,1,0);

        gl.glRotated(sim.getParameter().phi*180/Math.PI, 0,1,0);
        gl.glRotated(sim.getParameter().theta*180/Math.PI, 1,0,0);
        gl.glRotated(sim.getParameter().psi*180/Math.PI, 0,0,1);
        drawTop(gl);

        time += sim.getParameter().dt;
        sim.NDSolve(time, time + sim.getParameter().dt);
    }

    private void drawTop(GL2 gl){
        gl.glColor3d(0.8,0.2,0.2);
        glut.glutSolidCone(5,20,20,20);
    }

    private com.jogamp.opengl.util.gl2.GLUT glut = new com.jogamp.opengl.util.gl2.GLUT();

    @Override
    public void dispose(GLAutoDrawable drawable){ }
}
