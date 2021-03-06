/*
 * Extend particle to add things like the reference to the corresponding MC particle
 * Maybe also xF etc....
 * 
 */

package NovelFitters;

import org.jlab.clas.physics.Particle;


public class MyParticle extends Particle {	
	
	public int matchingMCPartIndex;
	public double m_chi2pid;
	
	//relatige to fastestd electron
	public float beta;
	public float FTOFPath;
	public float FTOFTime;
	public int FTOFsector;
	//0 means no PID
	public int PID;
	public MyParticle(int pid,float px,float py, float pz, float vx, float vy,float vz)
	{	
		super(pid,px,py,pz,vx,vy,vz);
		matchingMCPartIndex=-1;
	}	
}
