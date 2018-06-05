package NovelFitters;

import org.jlab.clas.physics.GenericKinematicFitter;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.PhysicsEvent;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataBank;
import org.jlab.clas.physics.*;

public class NovelBaseFitter extends GenericKinematicFitter {
	
	public Vector3 gNBoost;
	public double Walt;
	protected final Double m_beam;
	protected final boolean m_isMC;
	//indicate that the MC truth should be used
	protected final boolean m_useMC;
	protected double Q2;
	protected double x;
	protected double W;
	protected double nu;
	protected LorentzVector q;
	protected String bankName;
	
	
	
    public NovelBaseFitter(double beam, boolean isMC, boolean useMC) {
        super(beam);
        m_beam = beam;
        m_isMC=isMC;
        m_useMC=useMC;

    	if(useMC)
    		bankName=new String("MC::Particle");	
    	else
    		bankName=new String("REC::Particle");
        
    } 
    
    public double getQ2()
    {
    	return Q2;
    }
    public double getX()
    {
    	return x;
    }
    public double getW()
    {
    	return W;
    }
    public double getNu()
    {
    	return nu;
    }
    public LorentzVector getq()
    {
    	return q;
    }
    
    protected void computeKinematics(float px, float py, float pz)
    {
    	//
    		double m_e=0.511/1000.0;
    		double m_p=0.938;
    		LorentzVector lv_e=new LorentzVector();
    		LorentzVector lv_beam=new LorentzVector();
    		LorentzVector lv_target=new LorentzVector();
    		
    		lv_e.setPxPyPzM(px, py, pz, m_e);
    		lv_beam.setPxPyPzM(0.0,0.0,10.6,m_e);
    		lv_target.setPxPyPzM(0.0, 0.0, 0.0, m_p);
    		
    		q= new LorentzVector(lv_beam);
    		q.sub(lv_e);
    	//	System.out.println("q x " + q.px()+ " y: "+q.py() + " pz: " +q.pz() + " e: "+ q.e() );
    		//q= new LorentzVector(lv_e);
    		//q.sub(lv_beam);
    		Q2=(-1)*q.mass2();	
    		//need to multiply lorentz vectors... doesn't seem to be implemented in the jlab clas
    		//since this is lv_target*q /m_p and the momentum of the target is zero, it is 
    		//just the product of the energy parts, divided by m_p:
    		nu=lv_target.e()*q.e()/m_p;
   // 		System.out.println("target x " + lv_target.px()+ " y: "+lv_target.py() + " pz: " +lv_target.pz() + " e: "+ lv_target.e() );
    		x=Q2/(2*m_p*nu);
    		W=Math.sqrt(Q2*(1-x)/x);
    		LorentzVector gN=new LorentzVector(q);
    		gN.add(lv_target);
    	//	System.out.println("gN x " + gN.px()+ " y: "+gN.py() + " pz: " +gN.pz() + " e: "+ gN.e() );
    		Walt=gN.mass();
    		gNBoost=gN.boostVector();
    		gNBoost.negative();
    	//	System.out.println("x: "+x +" Q2 " + Q2 + " nu: "+ nu +"W: " + W);
    }
    /**
     * Returns PhysicsEvent object with reconstructed particles.
     *
     * @param event - DataEvent object
     * @return PhysicsEvent : event containing particles.
     */
    @Override
    public PhysicsEvent getPhysicsEvent(DataEvent event) 
    		{    
     	
    		boolean banks_test = true; // check to see if the event has all of the banks present
    		if (!(event.hasBank(bankName))) {
    			banks_test = false;
    			System.out.println("couldn't find bank" + bankName);
    		} 
    		if (banks_test) {
    			
    			int numElectrons=0;
    			PhysicsEvent physEvent = new PhysicsEvent();
    			HipoDataBank eventBank = (HipoDataBank) event.getBank(bankName); // load particle bank
    			for (int current_part = 0; current_part < eventBank.rows(); current_part++) {
    				int pid = eventBank.getInt("pid", current_part);
    				if(!m_isMC)
    				{
    					int status = eventBank.getInt("status", current_part);
    					float chi2pid = eventBank.getFloat("chi2pid", current_part);
    					//there seems to be a tendency to construct too many neutrons.
    					//the x*100 field in the status gives the number of scintillator responses
    					//good neutrons seem to have 2 
    					//also, even if there are two of those, the good one seems to have chi2pid < 999.0
    					//also, 4000 is CT and 2000 FT
    					if(pid == 2112 && (status< 2200 || status > 2999 || chi2pid >=999.0))
    						continue; 		
    				}
    				if (pid!=0 ) {
    						float vx = eventBank.getFloat("vx",current_part);
    						float vy = eventBank.getFloat("vy",current_part);
    						float vz = eventBank.getFloat("vz",current_part);
    						float px = eventBank.getFloat("px", current_part);
    						float py = eventBank.getFloat("py", current_part);
    						float pz = eventBank.getFloat("pz", current_part);
    					//System.out.println("pid: "+ pid +" pz: "+ pz +" status " + status + " chi2pid: " + chi2pid); 
    						if(pid==11 && numElectrons==0)
    	    					{
    							
    							//looks like the first one has the higher momentum, so probably the scattered one
    							numElectrons++;
    							double mom=Math.sqrt(px*px+py*py+pz*pz);
    							//System.out.println("found electron num: " + numElectrons+" mom: "+mom);
    							computeKinematics(px,py,pz);
    	    					}
    	    					 						
    						MyParticle part = new MyParticle(pid,px,py,pz,vx,vy,vz);
    						
    						physEvent.addParticle(part);
    				}
    						
    			}
            //check MC banks for a match
            	if(m_isMC){
            			
            	
            }	
            	return physEvent;
    		}    		  
    		return new PhysicsEvent(this.m_beam);
    }
    
   
}
