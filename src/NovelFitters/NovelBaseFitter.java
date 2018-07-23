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
	// indicate that the MC truth should be used
	protected final boolean m_useMC;
	protected double Q2;
	protected double x;
	protected double W;
	protected double nu;
	
	protected int foundLambda;
	protected LorentzVector q;
	protected LorentzVector lv_beam;
	protected LorentzVector lv_e;
	protected String bankName;
	protected static int pionLundCode=211;
	protected static int maxArrSize=400;
	protected static int[] part_Cal_PCAL_sector=new int[maxArrSize];
	protected static int[] calPindex=new int[maxArrSize];
	protected static double[] part_Cal_PCAL_x=new double[maxArrSize];
	protected static double[] part_Cal_PCAL_y=new double[maxArrSize];
	protected static double part_Cal_PCAL_E[]=new double[maxArrSize];
	
	protected static double part_Cal_CalTotal_E[]=new double[maxArrSize];
	
	protected static int[] part_DC_sector=new int[maxArrSize];
	protected static double[] part_DC_c1x=new double[maxArrSize];
	protected static double[] part_DC_c1y=new double[maxArrSize];

	protected static double[] part_DC_c2x=new double[maxArrSize];
	protected static double[] part_DC_c2y=new double[maxArrSize];
	protected static double[] part_DC_c3x=new double[maxArrSize];
	protected static double[] part_DC_c3y=new double[maxArrSize];
	
	protected static int[] FTOFHit=new int[maxArrSize];
	
	//local copy of the event, so we can use it in multiple methods
	protected DataEvent m_event;
	

	public NovelBaseFitter(double beam, boolean isMC, boolean useMC) {
		super(beam);
		m_beam = beam;
		m_isMC = isMC;
		m_useMC = useMC;

		if (useMC)
			bankName = new String("MC::Particle");
		else
			bankName = new String("REC::Particle");

	}

	public int getNumLambda()
	{
		return foundLambda;	
	}
	
	public double getQ2() {
		return Q2;
	}

	public double getX() {
		return x;
	}

	public double getW() {
		return W;
	}

	public double getNu() {
		return nu;
	}

	public LorentzVector getq() {
		return q;
	}

	public LorentzVector getL() {
		return lv_beam;
	}
	

	protected void computeKinematics(float px, float py, float pz) {
		//`
		double m_e = LundPID.Electron.mass();
		double m_p = LundPID.Proton.mass();
		lv_e = new LorentzVector();
		lv_beam = new LorentzVector();
		LorentzVector lv_target = new LorentzVector();

		lv_e.setPxPyPzM(px, py, pz, m_e);
		lv_beam.setPxPyPzM(0.0, 0.0, 10.6, m_e);
		lv_target.setPxPyPzM(0.0, 0.0, 0.0, m_p);

		q = new LorentzVector(lv_beam);
		q.sub(lv_e);
		// System.out.println("q x " + q.px()+ " y: "+q.py() + " pz: " +q.pz() + " e: "+
		// q.e() );
		// q= new LorentzVector(lv_e);
		// q.sub(lv_beam);
		Q2 = (-1) * q.mass2(); 
		// need to multiply lorentz vectors... doesn't seem to be implemented in the
		// jlab clas
		// since this is lv_target*q /m_p and the momentum of the target is zero, it is
		// just the product of the energy parts, divided by m_p:
		nu = lv_target.e() * q.e() / m_p;
		// System.out.println("target x " + lv_target.px()+ " y: "+lv_target.py() + "
		// pz: " +lv_target.pz() + " e: "+ lv_target.e() );
		x = Q2 / (2 * m_p * nu);
		W = Math.sqrt(Q2 * (1 - x) / x);
		LorentzVector gN = new LorentzVector(q);
		gN.add(lv_target);
		// System.out.println("gN x " + gN.px()+ " y: "+gN.py() + " pz: " +gN.pz() + "
		// e: "+ gN.e() );
		Walt = gN.mass();
		gNBoost = gN.boostVector();
		gNBoost.negative();
		
		LorentzVector lvQTst=new LorentzVector(q);
		LorentzVector lvPTst=new LorentzVector(lv_target);
		lvQTst.boost(gNBoost);
		lvPTst.boost(gNBoost);
		//System.out.println("qx: "+ lvQTst.px()+ " qy: " + lvQTst.py() + " qz: "+ lvQTst.pz() + " e: " + lvQTst.e() );
		//System.out.println("px: "+ lvPTst.px()+ " py: " + lvPTst.py() + " pz: "+ lvPTst.pz() + " e: " + lvPTst.e() );
		//System.out.println("px: "+ gNBoost.x()+ " py: " + gNBoost.y() + " pz: "+ gNBoost.z()  );
		// System.out.println("x: "+x +" Q2 " + Q2 + " nu: "+ nu +"W: " + W);
	}

	/**
	 * Returns PhysicsEvent object with reconstructed particles.
	 *
	 * @param event
	 *            - DataEvent object
	 * @return PhysicsEvent : event containing particles.
	 */
	@Override
	public PhysicsEvent getPhysicsEvent(DataEvent event) {
		foundLambda=0;	
		boolean banks_test = true; // check to see if the event has all of the banks present
		if (!(event.hasBank(bankName))) {
			banks_test = false;
			System.out.println("couldn't find bank" + bankName);
		}
		else
		{
			//System.out.println("bank_test fine");
		}
		m_event=event;
		if (banks_test) {

			int numElectrons = 0;
			PhysicsEvent physEvent = new PhysicsEvent();
			HipoDataBank eventBank = (HipoDataBank) event.getBank(bankName); // load particle bank
			for (int current_part = 0; current_part < eventBank.rows(); current_part++) {
				//System.out.println("get pid");
				int pid = eventBank.getInt("pid", current_part);
				//System.out.println("get status");
				int status=-1;
				float chi2pid=-1;
				try {
					//seems like mc does not have status or chi2pid
					if(!m_isMC)
					{
				 status = eventBank.getInt("status", current_part);
				
					chi2pid = eventBank.getFloat("chi2pid", current_part);
					}
				}
				catch(NullPointerException ex)
				{
					
				
				}
			
				//System.out.println("done");
				if (!m_isMC) {
					// there seems to be a tendency to construct too many neutrons.
					// the x*100 field in the status gives the number of scintillator responses
					// good neutrons seem to have 2
					// also, even if there are two of those, the good one seems to have chi2pid <
					// 999.0
					// also, 4000 is CT and 2000 FT
				if (pid == 2112 && (status < 2200 || status > 2999 || chi2pid >= 999.0))
						continue;
				}
				
				if (pid != 0) {
					float vx = eventBank.getFloat("vx", current_part);
					float vy = eventBank.getFloat("vy", current_part);
					float vz = eventBank.getFloat("vz", current_part);
					float px = eventBank.getFloat("px", current_part);
					float py = eventBank.getFloat("py", current_part);
					float pz = eventBank.getFloat("pz", current_part);
					// System.out.println("pid: "+ pid +" pz: "+ pz +" status " + status + "
					// chi2pid: " + chi2pid);
					if (pid == 11 && numElectrons == 0) {
						if(!survivesStefanElectronCuts()) {
							continue;	

						}
						
						
						
						
						// looks like the first one has the higher momentum, so probably the scattered
						// one
						numElectrons++;
						double mom = Math.sqrt(px * px + py * py + pz * pz);
						// System.out.println("found electron num: " + numElectrons+" mom: "+mom);
						computeKinematics(px, py, pz);
					}

					MyParticle part = new MyParticle(pid, px, py, pz, vx, vy, vz);
					part.m_chi2pid=chi2pid;
					physEvent.addParticle(part);
				}

			}
			// check MC banks for a match
			if (m_isMC) {
				if(event.hasBank("MC::Lund"))
				{
					//System.out.println("found lund bank");
					HipoDataBank eventBankLund = (HipoDataBank) event.getBank("MC::Lund"); // load particle bank
					for (int current_part = 0; current_part < eventBankLund.rows(); current_part++) {
						int pid = eventBankLund.getInt("pid", current_part);
						if(pid==3122)
						{
							foundLambda++;
							//System.out.println("Found lambda "+foundLambda);
						}
					}
				}

			}
			return physEvent;
		}
		return new PhysicsEvent(this.m_beam);
	}
	
	boolean survivesStefanElectronCuts()
	{
		
		//REC:particle
		//REC:calorimeter
		//REC:cherenkov
		//REC:scintillator
		//REC:trajectory
		
		if (!(event.hasBank(bankName))) {
			HipoDataBank eventBank = (HipoDataBank) event.getBank(bankName); // load particle bank
			for (int current_part = 0; current_part < eventBank.rows(); current_part++) {
				//System.out.println("get pid");
				int pid = eventBank.getInt("pid", current_part);
				//System.out.println("get status");
				int status=-1;
				float chi2pid=-1;
				try {
					//seems like mc does not have status or chi2pid
					if(!m_isMC)
					{
				 status = eventBank.getInt("status", current_part);	
			
		}
		
		
		
		return true;
	}
	
	int stefanHadronPID()
	{
	
		
		return pionLundCode;
	}
	
	//C++ code from stefan's note
	
	boolean EC_hit_position_fiducial_cut(int j){
	
		int sec_PCAL = part_Cal_PCAL_sector[j]-1;
		double x_PCAL = part_Cal_PCAL_x[j];
		double y_PCAL = part_Cal_PCAL_y[j];
		
		//pcal energy cut. electron leaves more 0.06 GeV in PCAL
		if(part_Cal_PCAL_E[j]<0.06)
			return false;
		
		//not sure if that is correct. I guess six setors * 6=2*180
		double Pival=0.5;
		double x_PCAL_rot = y_PCAL * Math.sin(sec_PCAL*60.0*Pival/180)
		+ x_PCAL * Math.cos(sec_PCAL*60.0*Pival/180);
		double y_PCAL_rot = y_PCAL * Math.cos(sec_PCAL*60.0*Pival/180)
		- x_PCAL * Math.sin(sec_PCAL*60.0*Pival/180);
		double angle_PCAL = 60;
		double height_PCAL = 45;
		double slope_PCAL = 1/Math.tan(0.5*angle_PCAL*Pival/180);
		double left_PCAL = (height_PCAL - slope_PCAL * y_PCAL_rot);
		double right_PCAL = (height_PCAL + slope_PCAL * y_PCAL_rot);
		double radius2_PCAL = Math.pow(height_PCAL+6,2)-Math.pow(y_PCAL_rot,2);
		if(x_PCAL_rot > left_PCAL && x_PCAL_rot > right_PCAL && Math.pow(x_PCAL_rot,2) > radius2_PCAL
		&& x_PCAL_rot < 372) return true;
		else return false;
		}
	boolean EC_sampling_fraction_cut(int j){
		double sigma_range = 3;
		double p0mean[] = {0.101621, 0.115216, 0.109593, 0.114007, 0.114176, 0.108219};
		double p1mean[] = {-0.1569, 0.129423, -0.200557, 0.00208234, -0.0478747, -0.236188};
		double p2mean[] = {0.0111515, 0.00903229, 0.0136201, 0.00959565, 0.0119958, 0.0191656};
		double p3mean[] = {-0.000966913, -0.000921184, -0.00130608, -0.000631906, -0.000827219, -0.00193565};
		double p0sigma[] = {0.0132654, 0.014592, 0.0211849, 0.0198346, 0.0176063, 0.0213921};
		double p1sigma[] = {0.00384941, 0.00340508, -0.0015, -0.000342043, 0.00119106, -0.00108573};
		double mean, double sigma, upper_lim_total , double lower_lim_total ;
		for(Int_t k = 0; k < 6; k++){
		if(part_Cal_PCAL_sector[j]-1 == k){
		mean = p0mean[k] *( 1 + part_p[j]/sqrt(pow(part_p[j],2) + p1mean[k])) + p2mean[k] * part_p[j] + p3mean[k] * pow(part_p[j],2);
		sigma = p0sigma[k] + p1sigma[k] * sqrt(part_p[j]);
		upper_lim_total = mean + sigma_range * sigma;
		lower_lim_total = mean - sigma_range * sigma;
		}
		}
		if(part_Cal_energy_total[j]/part_p[j] <= upper_lim_total && part_Cal_energy_total[j]/part_p[j] >= lower_lim_total) return true;
		else return false;
		}
	
	boolean DC_hit_position_region1_fiducial_cut(int j){
		double angle = 60;
		double height = 31;
		double Pival=0.5;
		int sec = part_DC_sector[j]-1;
		double x1_rot = part_DC_c1y[j] * sin(sec*60.0*Pival/180) + part_DC_c1x[j] * cos(sec*60.0*Pival/180);
		double y1_rot = part_DC_c1y[j] * cos(sec*60.0*Pival/180) - part_DC_c1x[j] * sin(sec*60.0*Pival/180);
		double slope = 1/tan(0.5*angle*Pival/180);
		double left = (height - slope * y1_rot);
		double right = (height + slope * y1_rot);
		double radius2_DCr1 = Math.pow(32,2)-Math.pow(y1_rot,2);
		if (x1_rot > left && x1_rot > right && Math.pow(x1_rot,2) > radius2_DCr1) return true;
		else return false;
		}
	
	
	boolean loadFTOF()
	{
		String bankName="REC::Scintillator";
		if (!(m_event.hasBank(bankName))) {
			System.out.println("couldn't find bank" + bankName);
			return false;
		}
		else
		{
			HipoDataBank eventBank = (HipoDataBank) m_event.getBank(bankName);
			int counter=0;
			for (int current_part = 0; current_part < eventBank.rows(); current_part++) {
				//System.out.println("get pid");
				int pindex = eventBank.getInt("pindex", current_part);
				int detID = eventBank.getInt("detector", current_part);
				if(pindex>=maxArrSize)
				{
					System.out.println("too many particles is trajectory bank: "+pindex);
					break;
				}
				//entrance point to region 1
				if(detID==12)
				{
					FTOFHit[pindex]=1;
				}
				else
				{
					FTOFHit[pindex]=0;
				}
			
			}
		}
		return false;
	
	}
	
	
	boolean loadDCInfo()
	{
		String bankName="REC::Trajectory";
		if (!(m_event.hasBank(bankName))) {
			System.out.println("couldn't find bank" + bankName);
			return false;
		}
		else
		{
			HipoDataBank eventBank = (HipoDataBank) m_event.getBank(bankName);
			int counter=0;
			for (int current_part = 0; current_part < eventBank.rows(); current_part++) {
				//System.out.println("get pid");
				int pindex = eventBank.getInt("pindex", current_part);
				int detID = eventBank.getInt("detID", current_part);
				if(pindex>=maxArrSize)
				{
					System.out.println("too many particles is trajectory bank: "+pindex);
					break;
				}
				//entrance point to region 1
				if(detID==12)
				{
					int sector = eventBank.getInt("sector", current_part);	
					part_DC_sector[pindex]=sector;
					double x = eventBank.getDouble("x", current_part);
					part_DC_c1x[pindex]=x;
					double y = eventBank.getDouble("y", current_part);
					part_DC_c1y[pindex]=y;
				}
				if(detID==24)
				{
					double x = eventBank.getDouble("x", current_part);
					part_DC_c2x[pindex]=x;
					double y = eventBank.getDouble("y", current_part);
					part_DC_c2y[pindex]=y;
				}
			if(detID==36)
			{
					double x = eventBank.getDouble("x", current_part);
					part_DC_c2x[pindex]=x;
					double y = eventBank.getDouble("y", current_part);
					part_DC_c2y[pindex]=y;
					
					counter++;
				}
			}
		}
		return false;
	
	}
	
	
	
	
	
	boolean loadCalInfo()
	{
		
		String bankName="REC::Calorimeter";
		if (!(m_event.hasBank(bankName))) {
			System.out.println("couldn't find bank" + bankName);
			return false;
		}
		else
		{
			
		for(int i=0;i<maxArrSize;i++)
		{
			
		}
			HipoDataBank eventBank = (HipoDataBank) m_event.getBank(bankName);
			int counter=0;
			for (int current_part = 0; current_part < eventBank.rows(); current_part++) {
				//System.out.println("get pid");
				int pindex = eventBank.getInt("pindex", current_part);
				if(pindex>=maxArrSize)
				{
					System.out.println("too many particles is calorimeter bank: "+pindex);
					break;
				}
				int layer=eventBank.getInt("layer", current_part);
				//pcal is layer 1, 4 ECinner, 7 ECOuter
				double energy=eventBank.getDouble("energy",current_part);
				if(layer==1)
				{
					calPindex[counter]=pindex;
					int sector = eventBank.getInt("sector", current_part);
					part_Cal_PCAL_sector[pindex]=sector;
					double x = eventBank.getDouble("x", current_part);
					part_Cal_PCAL_x[pindex]=x;
					double y = eventBank.getDouble("y", current_part);
					part_Cal_PCAL_y[pindex]=y;
					part_Cal_PCAL_E[pindex]=energy;		
				}
				if(layer==1 || layer==4 || layer==7)
				{
					part_CalTotal+=energy
				}
				counter++;
			}
			
		}
		return false;
	}
		

	
			boolean hasPCALInfo(int particleIndex)
			{
				string bankName="REC::Calorimeter";
			
				if (!(m_event.hasBank(bankName))) {
					
					System.out.println("couldn't find bank" + bankName);
					return false;
				}
				else
				{
					HipoDataBank eventBank = (HipoDataBank) m_event.getBank(bankName);
					for (int current_part = 0; current_part < eventBank.rows(); current_part++) {
						//System.out.println("get pid");
						int pindex = eventBank.getInt("pindex", current_part);
						if(pindex==particleIndex)
							return true;
					}
				}
				return false;
			}
	
			
			

}
