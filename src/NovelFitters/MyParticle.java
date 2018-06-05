/*
 * Extend particle to add things like the reference to the corresponding MC particle
 * Maybe also xF etc....
 * 
 */

package NovelFitters;

import org.jlab.clas.physics.Particle;


public class MyParticle extends Particle {	
	
	public int matchingMCPartIndex;
	
	public MyParticle(int pid,float px,float py, float pz, float vx, float vy,float vz)
	{	
		super(pid,px,py,pz,vx,vy,vz);
		matchingMCPartIndex=0;
	}	
}
