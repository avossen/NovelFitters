package NovelFitters;

import org.jlab.clas.physics.GenericKinematicFitter;
import org.jlab.groot.data.H2F;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.PhysicsEvent;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataBank;

import java.util.Arrays;

import org.jlab.clas.physics.*;

import org.apache.commons.math3.*;

public class NovelBaseFitter extends GenericKinematicFitter {
public static int debugEvent=48478732;

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
	protected double y;
	protected float torus;
	protected float solenoid;
	protected int runNumber;
	protected int evtNumber;
	protected float electronTime;
	protected boolean population_weighting = true;
	protected boolean population_weighting_CD = true;
	protected boolean outbending = false;

	protected int helicity;
	protected int foundLambda;
	protected LorentzVector q;
	protected LorentzVector lv_beam;
	protected LorentzVector lv_e;
	protected String bankName;
	protected String bankNameMC;
	public static boolean useStefanPIDCuts=false;
	public static boolean useStefanElectronCuts = false;
	public static boolean noVertexCut=false;
	public static boolean forceEBElectronId=true;
	public static boolean useAllStefanElectronCuts=false;
	// for rec software before June, we need to access timebasedtracks bank
	// to get sector. After that it is in in REC::Tracks
	public static boolean useTimeBasedTracks = false;

	public static boolean useStefanHadronCuts = false;
	protected static int maxArrSize = 400;
	// don't make the arrays static, since we will have several instances of the
	// fitter running (e.g. data and mc)
	protected int[] part_Cal_PCAL_sector = new int[maxArrSize];
	protected int[] calPindex = new int[maxArrSize];
	protected float[] part_Cal_PCAL_x = new float[maxArrSize];
	protected float[] part_Cal_PCAL_y = new float[maxArrSize];
	protected float[] PCAL_lu = new float[maxArrSize];
	protected float[] PCAL_lv = new float[maxArrSize];
	protected float[] PCAL_lw = new float[maxArrSize];
	
	
	protected float part_Cal_PCAL_E[] = new float[maxArrSize];
	protected int trkSectors[] = new int[maxArrSize];
	protected double part_Cal_CalTotal_E[] = new double[maxArrSize];

	protected int[] part_DC_sector = new int[maxArrSize];
	protected float[] part_DC_c1x = new float[maxArrSize];
	protected float[] part_DC_c1y = new float[maxArrSize];

	protected float[] part_DC_c2x = new float[maxArrSize];
	protected float[] part_DC_c2y = new float[maxArrSize];
	protected float[] part_DC_c3x = new float[maxArrSize];
	protected float[] part_DC_c3y = new float[maxArrSize];
	protected float[] part_p = new float[maxArrSize];
	protected float[] part_beta = new float[maxArrSize];
	protected float[] part_vz = new float[maxArrSize];
	protected int[] part_charge = new int[maxArrSize];

	protected int[] FTOFHit = new int[maxArrSize];
	protected int[] FTOFSector = new int[maxArrSize];
	protected float[] FTOFTime = new float[maxArrSize];
	protected float[] FTOFPath = new float[maxArrSize];

	protected H2F[] dcAcceptedPos = new H2F[3];
	protected H2F[] dcAllPos = new H2F[3];

	// local copy of the event, so we can use it in multiple methods
	protected DataEvent m_event;
	protected MyParticle scatteredElectron;

	public NovelBaseFitter(double beam, boolean isMC, boolean useMC) {
		super(beam);
		m_beam = beam;
		m_isMC = isMC;
		m_useMC = useMC;

		for (int i = 0; i < 3; i++) {
			String title = "dc" + (i + 1) + "AllPos";
			dcAllPos[i] = new H2F(title, title, 100, -300, 300, 100, -300, 300);
			title = "dc" + (i + 1) + "AcceptedPos";
			dcAcceptedPos[i] = new H2F(title, title, 100, -300, 300, 100, -300, 300);
		}

		bankNameMC = new String("MC::Particle");
		bankName = new String("REC::Particle");
		// if (useMC)
		// bankName = new String("MC::Particle");
		// else
		// bankName = new String("REC::Particle");

	}

	public int getNumLambda() {
		return foundLambda;
	}

	public int getBeamHelicity() {
		return helicity;
	}

	public float getTorus()
	{
		return torus;
	}
	public float getSolenoid()
	{
		return solenoid;
	}
	
	public double getQ2() {
		return Q2;
	}
	
	public int getRunNumber()
	{
		return runNumber;
	}

	public int getEvtNumber()
	{
		return evtNumber;
	}
	public double getY() {
		return y;
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

	protected void cleanArrays() {
		Arrays.fill(this.trkSectors, 0);
		Arrays.fill(part_Cal_PCAL_sector, 0);
		Arrays.fill(calPindex, 0);
		Arrays.fill(part_Cal_PCAL_x, (float) 0.0);
		Arrays.fill(part_Cal_PCAL_y, (float) 0.0);
		Arrays.fill(PCAL_lu, (float) 0.0);
		Arrays.fill(PCAL_lv, (float) 0.0);
		Arrays.fill(PCAL_lw, (float) 0.0);
		
		
		Arrays.fill(part_Cal_PCAL_E, (float) 0.0);

		Arrays.fill(part_Cal_CalTotal_E, 0.0);
		Arrays.fill(part_DC_sector, 0);
		Arrays.fill(part_DC_c1x, (float) 0.0);
		Arrays.fill(part_DC_c1y, (float) 0.0);

		Arrays.fill(part_DC_c2x, (float) 0.0);
		Arrays.fill(part_DC_c2y, (float) 0.0);

		Arrays.fill(part_DC_c3x, (float) 0.0);
		Arrays.fill(part_DC_c3y, (float) 0.0);
		Arrays.fill(part_p, (float) 0.0);
		Arrays.fill(part_beta, (float) 0.0);
		Arrays.fill(part_vz, (float) 0.0);
		Arrays.fill(part_charge, 0);
		Arrays.fill(FTOFHit, 0);
		Arrays.fill(FTOFTime, (float) -1.0);
		Arrays.fill(FTOFPath, (float) -1.0);
		Arrays.fill(FTOFSector, -1);
		helicity = 7;
	}

	protected void computeKinematics(float px, float py, float pz) {
		// `
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
		y = q.e() / lv_beam.e();
		// System.out.println("target x " + lv_target.px()+ " y: "+lv_target.py() + "
		// pz: " +lv_target.pz() + " e: "+ lv_target.e() );
		x = Q2 / (2 * m_p * nu);
		W = Math.sqrt(m_p+Q2 * (1 - x) / x);
		LorentzVector gN = new LorentzVector(q);
		gN.add(lv_target);
		// System.out.println("gN x " + gN.px()+ " y: "+gN.py() + " pz: " +gN.pz() + "
		// e: "+ gN.e() );
		Walt = gN.mass();
		//System.out.println("W: " + W  + " Walt: "+ Walt);
		
		gNBoost = gN.boostVector();
		gNBoost.negative();

		LorentzVector lvQTst = new LorentzVector(q);
		LorentzVector lvPTst = new LorentzVector(lv_target);
		lvQTst.boost(gNBoost);
		lvPTst.boost(gNBoost);
		// System.out.println("qx: "+ lvQTst.px()+ " qy: " + lvQTst.py() + " qz: "+
		// lvQTst.pz() + " e: " + lvQTst.e() );
		// System.out.println("px: "+ lvPTst.px()+ " py: " + lvPTst.py() + " pz: "+
		// lvPTst.pz() + " e: " + lvPTst.e() );
		// System.out.println("px: "+ gNBoost.x()+ " py: " + gNBoost.y() + " pz: "+
		// gNBoost.z() );
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
		
		this.cleanArrays();
		foundLambda = 0;
		boolean banks_test = true; // check to see if the event has all of the banks present

		//if we use the mc truth, we don't  do fiducial cuts etc...
		if(m_useMC)
		{
			useStefanElectronCuts=false;
			useStefanHadronCuts=false;
			useStefanPIDCuts=false;
		}
		
		
		if (!(event.hasBank(bankName)) && (!m_useMC || event.hasBank(bankNameMC))) {
			banks_test = false;
			// System.out.println("couldn't find bank" + bankName);
		} else {
			// System.out.println("bank_test fine");
		}
		m_event = event;

		boolean loadedDC = this.loadDCInfo();

		boolean loadedFTOF = this.loadFTOF();
		boolean loadedCal = this.loadCalInfo();

		try {
			//we only care about loading the detector info if we use stefan's cuts
			//otherwise all we need should be in the particle bank
			if ((!banks_test) || ((useStefanElectronCuts || useStefanHadronCuts) &&(!loadedDC || !loadedFTOF || !loadedCal))) {
				throw new Exception("bank missing");
			}
			HipoDataBank eventBank = (HipoDataBank) event.getBank(bankName); // load particle bank
			HipoDataBank eventBankMC = null;

			if (m_useMC)
			{
				eventBankMC = (HipoDataBank) event.getBank(bankNameMC);
				//eventBank now points to evenBankMC
				eventBank=eventBankMC;
			}

			try {
				
				helicity = getHelicityAndSetRunEvtNumbers();
			
				if(this.evtNumber==debugEvent)
				{
					System.out.println("found our run.. (novel base)");
				}
			//	System.out.println("helicity: " + helicity);

				// seems like mc does not have status or chi2pid
				
			} catch (NullPointerException ex) {
			}
			// System.out.println("check for electron");
			if (!findScatteredElectron(eventBank)) {
				if(this.evtNumber==debugEvent)
				{
					System.out.println("no electron.. (novel base)");
				}
				throw new Exception("no electron");
			}
			// System.out.println("got electron");
			PhysicsEvent physEvent = new PhysicsEvent();
			//should have been set already by the findScatteredElectron routine
			//this.scatteredElectron.PID=11;
			physEvent.addParticle(this.scatteredElectron);
			// if(eventBank.rows()>=maxArrSize)
			// continue;

			// go for the hadrons-->since in the mc case this points now to the mc
			//bank, we would run over all mc particles...
			for (int current_part = 0; current_part < eventBank.rows(); current_part++) {
				// if(eventBank.rows()>5)
				// System.out.println("num particles: " +eventBank.rows());
				if (current_part >= maxArrSize) {
					System.out.println("max array size..");
					break;
				}
				// System.out.println("get pid");
				int pid = eventBank.getInt("pid", current_part);

				// System.out.println("get status");
				int status = -1;
				float chi2pid = -1;

				try {
				
					if (!m_isMC ) {
						status = eventBank.getInt("status", current_part);
						chi2pid = eventBank.getFloat("chi2pid", current_part);
					}
				} catch (NullPointerException ex) {
				}

				// System.out.println("done");
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

				//so either we get event builder PID or we use stefan's PID
				if (pid != 0 || (useStefanHadronCuts && useStefanPIDCuts) ) {
					// for now only interested in pions and kaons. But use all particles in the MC
					// that means that a particle identified as a pion might loose its MC partner if
					// that is not a pion
					
					//not the most elegant, but we need a pid!=0 to construct a particle
					//that we can pass to the PID stuff, so we just assume a pion here
					//this is only possible if we use Stefan's cuts for hadron and pid,
					//so for that we reassess the pid anyways
					if(pid==0)
						pid=211;
					
					
					float vx = eventBank.getFloat("vx", current_part);
					float vy = eventBank.getFloat("vy", current_part);
					float vz = eventBank.getFloat("vz", current_part);
					float px = eventBank.getFloat("px", current_part);
					float py = eventBank.getFloat("py", current_part);
					float pz = eventBank.getFloat("pz", current_part);
					float beta=0;
					if(!m_useMC)
					 beta = eventBank.getFloat("beta", current_part);

					

					// System.out.println("pid: "+ pid +" pz: "+ pz +" status " + status + "
					// chi2pid: " + chi2pid);

					// vertex cuts for hadrons? probably not...
					// if(vz< -7 || vz>10)
					// continue;

					float mom = (float) Math.sqrt(px * px + py * py + pz * pz);
					// the pindex of in the particle table should just be the index
					this.part_p[current_part] = mom;
					this.part_beta[current_part] = beta;
					part_vz[current_part] = vz;

					MyParticle part = null;
					// the idea of the below is that we only add particles to the MC events that
					// are accepted in the reconstruction. That way we can do the matching but don't
					// store pairs that are unnceessary
					// this would be done differently if we want to reconsruct 'true' MC asymmetries
				
					//since eventBank and eventBankMC is  now the same, don't need the test anymore
						part = new MyParticle(pid, px, py, pz, vx, vy, vz);
				
					part_charge[current_part] = part.charge();
					//initialize for MC case
					part.FTOFsector=-1;
					part.FTOFTime=-1;
					part.FTOFPath=-1;
					part.m_chi2pid=-1;
					if(!m_useMC)
					{
						part.FTOFsector = this.FTOFSector[current_part];
						// System.out.println("ftofsector: " + FTOFSector[current_part]);
						if (this.FTOFSector[current_part] >= 0) {
							part.FTOFTime = this.FTOFTime[current_part] - this.electronTime;
							part.FTOFPath = this.FTOFPath[current_part];
							float speedOLight = (float) 29.9792; // in cm per ns as per units in hip bank
							// according to Stefan's code, this should be 10^7, but that makes beta too
							// large, maybe speed of light is not in the right units
							// part.beta=(float)Math.pow(10,1)*part.FTOFPath/(part.FTOFTime*speedOLight);
							part.beta = beta;
							// System.out.println("beta: " + beta);
						/// this is the computed one, but there seems to be somthing
							// wrongpart.FTOFPath/(part.FTOFTime*speedOLight);
							// System.out.println("beta computed: " + part.beta + " beta from bank " +
							// beta);
						} else {
							part.FTOFTime = -1;
							part.beta = -1;
						}
						part.m_chi2pid = chi2pid;
					}
					int myPid = pid;
					//these are just the fiducial cuts
					if(useStefanHadronCuts)
					{
						if(0==stefanHadronFiducialCuts(current_part,part))
						{
							continue;
						}
					}	
					if (useStefanHadronCuts && useStefanPIDCuts) {
						myPid = this.stefanHadronPID(current_part, part);
						
					}
					if (myPid == 0)
						continue;
					if(this.evtNumber==debugEvent && Math.abs((myPid))==211)
					{
						System.out.println("found pion ("+myPid+") with mom x: " + px + " py " + py+ " pz "+ pz);

					}
					if(!m_useMC)
					{
						if (Math.abs(pid) != LundPID.Pion.lundCode() && Math.abs(pid) != LundPID.Kaon.lundCode())
							continue;
					}
					else
					{
						//for MC we should basically take everything to not 
						//bias the matching. Exception are neutrons which can be plentiful
						if(Math.abs(pid)==2112)
							continue;
					}
									
					// System.out.println("pid: " + pid + " stefan's id "+myPid);
					// System.out.println("Adding particle");
					if(this.evtNumber==debugEvent)
					{
						System.out.println("adding "+ pid);
					}
					physEvent.addParticle(part);
				}

			}

			// check MC banks for a match
			if (m_isMC) {
				if (event.hasBank("MC::Lund")) {
					// System.out.println("found lund bank");
					HipoDataBank eventBankLund = (HipoDataBank) event.getBank("MC::Lund"); // load particle bank
					for (int current_part = 0; current_part < eventBankLund.rows(); current_part++) {
						int pid = eventBankLund.getInt("pid", current_part);
						if (pid == 3122) {
							foundLambda++;
							// System.out.println("Found lambda "+foundLambda);
						}
					}
				}

			}
			// System.out.println("returning phys event with " + physEvent.count() + "
			// particles");
			return physEvent;
		} catch (Exception e) {
			// System.out.println("exception: "+ e.getMessage());
			return new PhysicsEvent(this.m_beam);
		}

	}

	boolean findScatteredElectron(HipoDataBank eventBank) {
		boolean foundElectron = false;
		float maxMom = -1000;
		float maxPx = 0;
		float maxPy = 0;
		float maxPz = 0;

		float maxVx = 0;
		float maxVy = 0;
		float maxVz = 0;

		for (int current_part = 0; current_part < eventBank.rows(); current_part++) {

			int pid = eventBank.getInt("pid", current_part);
			// System.out.println("get status");
			if(forceEBElectronId&&pid!=LundPID.Electron.lundCode())
				continue;
			if (pid != 0) {
				if(this.evtNumber==debugEvent)
				{
					System.out.println("looking for e, found pid " + pid);
				}
				float vx = eventBank.getFloat("vx", current_part);
				float vy = eventBank.getFloat("vy", current_part);
				float vz = eventBank.getFloat("vz", current_part);
				float px = eventBank.getFloat("px", current_part);
				float py = eventBank.getFloat("py", current_part);
				float pz = eventBank.getFloat("pz", current_part);
				float chi2pid=-1;
				try {
					if(!m_isMC)
					{
						chi2pid=eventBank.getFloat("chi2pid", current_part);
					}
			} catch (NullPointerException ex) {
			}
				if(Math.abs(chi2pid)>=9999.0)
					continue;
				// System.out.println("pid: "+ pid +" pz: "+ pz +" status " + status + "
				// chi2pid: " + chi2pid);
				if(this.evtNumber==debugEvent)
				{
					System.out.println("testing vz: " + vz);
				}
				/// Stefan's vz cut is sector dependent. Here only rough
				if (useStefanElectronCuts&& !DC_z_vertex_cut(current_part, vz))
				{
					if(this.evtNumber==debugEvent)
					{
						System.out.println("vtx cut");
					}
					continue;
				}
				// System.out.println("after vertex cut");

				float mom = (float) Math.sqrt(px * px + py * py + pz * pz);
				// the pindex of in the particle table should just be the index
				this.part_p[current_part] = mom;
				// System.out.println("electron mom: " + mom);
				// trigger thresholds is 1.5 GeV
				// if (pid == LundPID.Electron.lundCode() && mom>1.5) {
				
				//for x-check purposes
				if(mom< 2.0)
				{
					if(this.evtNumber==debugEvent)
					{
						System.out.println("mom cut " + mom);
					}
					continue;
				}
					//use the 1.5 cut also for the 'no stefan' case
				if(!m_useMC &&!useStefanElectronCuts && mom < 1.5)
					continue;
				
				if (useStefanElectronCuts &&mom < 1.5) 
					continue;
					// System.out.println("after mom cut");
					if (useStefanElectronCuts && !survivesStefanElectronCuts(current_part)) {
						continue;
					}
					else
					{
						
					}
					// if we don't use stefan's cut, use default assignment
					if (!useStefanElectronCuts && !(pid == LundPID.Electron.lundCode()))
						continue;
					// System.out.println("found electron");
					foundElectron = true;
					if (mom > maxMom) {
						maxMom = mom;
						maxPx = px;
						maxPy = py;
						maxPz = pz;
						if(this.evtNumber==debugEvent)
						{
							System.out.println("found electron with mom x: " + px + " py " + py+ " pz "+ pz);
	
						}
						maxVx = vx;
						maxVy = vy;
						maxVz = vz;
					}

					// take this electron to define FTOF start time
					this.electronTime = this.FTOFTime[current_part];
					// looks like the first one has the higher momentum, so probably the scattered
					// one

					// System.out.println("found electron with px " + px + " py "+ py + " pz " +
					// pz);
					// System.out.println("found electron num: " + numElectrons+" mom: "+mom);

				
			}
		}
		if (foundElectron) {
			int pid = LundPID.Electron.lundCode();
			this.scatteredElectron = new MyParticle(pid, maxPx, maxPy, maxPz, maxVx, maxVy, maxVz);
			computeKinematics(maxPx, maxPy, maxPz);
		}
		return foundElectron;
	}

	boolean survivesStefanElectronCuts(int partIndex) {

		// REC:particle
		// REC:calorimeter
		// REC:cherenkov
		// REC:scintillator
		// REC:trajectory

		boolean ecFiducialCuts = EC_hit_position_fiducial_cut(partIndex);
		//boolean ecFiducialCuts = EC_hit_position_fiducial_cut_natural(partIndex);
		boolean dcFiducialCuts = DC_hit_position_fiducial_cut(partIndex, false, false);
		//boolean ecEnergyDeposit = EC_sampling_fraction_cut(partIndex);
		boolean ecEnergyDeposit=true;
		boolean hasFtofHit = false;
		if (this.FTOFHit[partIndex] > 0)
			hasFtofHit = true;

		// pcal energy
		boolean pcalECut = true;
		if (this.part_Cal_PCAL_E[partIndex] < 0.06)
			pcalECut = false;
		
		
		
		
		// if(dcFiducialCuts)
		// System.out.println("dc fid cuts!");
		// System.out.println("ec cut: " + ecFiducialCuts + " dc "+dcFiducialCuts + " ec
		// energy "+ ecEnergyDeposit + " has ftof "+hasFtofHit+" pcal "+ pcalECut);
		//kick out the ecal and ftof hits (so everything besides fiducial)
		
 		
		if(this.evtNumber==debugEvent)
		{
			System.out.println(debugEvent+" survives fid cuts? "+ ecFiducialCuts + " dc: "+ dcFiducialCuts);
		}
		if(useAllStefanElectronCuts)
			return (ecFiducialCuts && dcFiducialCuts && ecEnergyDeposit && hasFtofHit && pcalECut);
		else
			return (ecFiducialCuts && dcFiducialCuts );
	}

	void saveHitHistograms()
	{
		EmbeddedCanvas can = new EmbeddedCanvas();
		can.setSize(1200, 600);
		can.divide(3, 2);
		//can_piPi.setAxisTitleSize(24);
		//can_piPi.setAxisFontSize(24);

	//	can_piPi.setTitleSize(24);
		can.cd(0);
		can.getPad(0).getAxisX().setTitle("dc1 x");
		can.getPad(0).getAxisY().setTitle("dc1 y");
		can.draw(this.dcAllPos[0]);
		can.cd(1);
		can.getPad(1).getAxisX().setTitle("dc2 x");
		can.getPad(1).getAxisY().setTitle("dc2 y");
		can.draw(this.dcAllPos[1]);
		can.cd(2);
		can.getPad(2).getAxisX().setTitle("dc3 x");
		can.getPad(2).getAxisY().setTitle("dc3 y");
		can.draw(this.dcAllPos[2]);
		
		can.cd(3);
		can.getPad(3).getAxisX().setTitle("dc1 x");
		can.getPad(3).getAxisY().setTitle("dc1 y");
		can.draw(this.dcAcceptedPos[0]);
		can.cd(4);
		can.getPad(4).getAxisX().setTitle("dc2 x");
		can.getPad(4).getAxisY().setTitle("dc2 y");
		can.draw(this.dcAcceptedPos[1]);
		can.cd(5);
		can.getPad(5).getAxisX().setTitle("dc3 x");
		can.getPad(5).getAxisY().setTitle("dc3 y");
		can.draw(this.dcAcceptedPos[2]);
		
		can.save("hitmap.png");

		
	}
	
	int stefanHadronFiducialCuts(int partIndex, MyParticle part)
	{
		boolean isPos = false;
		if (part.charge() > 0)
			isPos = true;
		boolean dcFiducialCut = DC_hit_position_fiducial_cut(partIndex, true, isPos);
		if (!dcFiducialCut)
			return 0;
		return 1;
	}
	int stefanHadronPID(int partIndex, MyParticle part) {
		boolean isPos = false;
		if (part.charge() > 0)
			isPos = true;
		
		int chargeFactor = (-1);
		if (isPos)
			chargeFactor = 1;

		// return (chargeFactor)*LundPID.Pion.lundCode();
		if (maximum_probability_cut(partIndex, (chargeFactor) * LundPID.Pion.lundCode(), 0.27, 99.73)) {
			if ((isPos && pip_delta_vz_cut(partIndex)) || (!isPos && pim_delta_vz_cut(partIndex)))
				return (chargeFactor) * LundPID.Pion.lundCode();
		}
		if (maximum_probability_cut(partIndex, (chargeFactor) * LundPID.Kaon.lundCode(), 0.27, 99.73)) {
			if ((isPos && Kp_delta_vz_cut(partIndex)) || (!isPos && Km_delta_vz_cut(partIndex)))
				return (chargeFactor) * LundPID.Kaon.lundCode();
		}
		if (maximum_probability_cut(partIndex, (chargeFactor) * LundPID.Proton.lundCode(), 0.27, 99.73)) {
			if ((isPos && prot_delta_vz_cut(partIndex)) || (!isPos && prot_delta_vz_cut(partIndex)))
				return (chargeFactor) * LundPID.Proton.lundCode();
		}

		// nothing identified
		return 0;
	}
//for electron
	boolean DC_z_vertex_cut(int j, float vz) {
if(noVertexCut)
	return true;
		double vz_min_sect[] = {-12, -12, -12, -14, -12, -12};
		double vz_max_sect[] = {10, 10, 10, 8, 10, 10};	
		double vz_min = 0;
		double vz_max = 0;

		for (int k = 0; k < 6; k++) {
			// uses pcal... I guess fine for electrons where we demand a pcal hit anyways.
			if (part_Cal_PCAL_sector[j] - 1 == k) {
				vz_min = vz_min_sect[k];
				vz_max = vz_max_sect[k];
			}
		}

		if (vz > vz_min && vz < vz_max)
			return true;
		else
			return false;
	}

	
	// C++ code from stefan's note
	
	
	//cut in the natural pcal co-ordinates-->outdated numbers
	boolean EC_hit_position_fiducial_cut_natural(int j){
		double u = PCAL_lu[j];
		double v = PCAL_lv[j];
		double w = PCAL_lw[j];
		double min_u = 8;
		double max_u = 400;
		double min_v = 8;
		double max_v = 400;
		double min_w = 8;
		double max_w = 400;
		if(u > min_u && u < max_u && v > min_v && v < max_v && w > min_w && w < max_w) return true;
		else return false;
		}
	boolean EC_hit_position_fiducial_cut(int j) {

		int sec_PCAL = part_Cal_PCAL_sector[j] - 1;
		double x_PCAL = part_Cal_PCAL_x[j];
		double y_PCAL = part_Cal_PCAL_y[j];

		// pcal energy cut. electron leaves more 0.06 GeV in PCAL
		//we test this extra, this function is only fiducial cuts
	//	if (part_Cal_PCAL_E[j] < 0.06)
	//		return false;

		// not sure if that is correct. I guess six setors * 6=2*180
		double Pival = Math.PI;
		double x_PCAL_rot = y_PCAL * Math.sin(sec_PCAL * 60.0 * Pival / 180)
				+ x_PCAL * Math.cos(sec_PCAL * 60.0 * Pival / 180);
		double y_PCAL_rot = y_PCAL * Math.cos(sec_PCAL * 60.0 * Pival / 180)
				- x_PCAL * Math.sin(sec_PCAL * 60.0 * Pival / 180);
		double angle_PCAL = 60;
		double height_PCAL = 47;
		double slope_PCAL = 1 / Math.tan(0.5 * angle_PCAL * Pival / 180);
		double left_PCAL = (height_PCAL - slope_PCAL * y_PCAL_rot);
		double right_PCAL = (height_PCAL + slope_PCAL * y_PCAL_rot);
		double radius2_PCAL = Math.pow(height_PCAL + 6, 2) - Math.pow(y_PCAL_rot, 2);
		if (x_PCAL_rot > left_PCAL && x_PCAL_rot > right_PCAL && Math.pow(x_PCAL_rot, 2) > radius2_PCAL
				&& x_PCAL_rot < 372)
			return true;
		else
			return false;
	}

	boolean EC_sampling_fraction_cut(int j) {
		double sigma_range = 3;
		double p0mean[] = {0.106333, 0.113711, 0.107714, 0.113276, 0.115548, 0.11108};
		double p1mean[] = {-0.374003, 0.164037, -0.101566, 0.22524,0.272903, 0.0370852};
		double p2mean[] = {0.00816235, 0.00390166, 0.00832663, 0.00324039,
				0.00376747, 0.00899919};
		double p3mean[] = {-0.000789648, -0.000432195, -0.00091734,
				-0.00013304, -0.000173935, -0.000932962};
		double p0sigma[] = {0.0162041, 0.0256151, 0.00996036, 0.0174414,
				0.0195056, 0.0115662};
		double p1sigma[] = {0.00472716, -0.00465669, 0.0140508, 0.00455405,
				0.000429308, 0.0119683};
		double mean = 0.0;
		double sigma = 0.0;
		double upper_lim_total = 0.0;
		double lower_lim_total = 0.0;
		for (int k = 0; k < 6; k++) {
			if (part_Cal_PCAL_sector[j] - 1 == k) {
				mean = p0mean[k] * (1 + part_p[j] / Math.sqrt(Math.pow(part_p[j], 2) + p1mean[k]))
						+ p2mean[k] * part_p[j] + p3mean[k] * Math.pow(part_p[j], 2);
				sigma = p0sigma[k] + p1sigma[k] * Math.sqrt(part_p[j]);
				upper_lim_total = mean + sigma_range * sigma;
				lower_lim_total = mean - sigma_range * sigma;
			}
		}
		if (this.part_Cal_CalTotal_E[j] / part_p[j] <= upper_lim_total
				&& this.part_Cal_CalTotal_E[j] / part_p[j] >= lower_lim_total)
			return true;
		else
			return false;
	}

	// these cuts are obviously dependendent on inbending/outbending
	boolean 	DC_hit_position_fiducial_cut(int j, boolean isHadron, boolean positive) {
		double Pival = Math.PI;
		double angle = 60;
		boolean cut[] = new boolean[3];

		double heightElectron[] = { 29, 43, 48 };
		double radiusElectron[] = { 30, 45, 50 };

		double heightHPlus[] = { 19, 38, 70 };
		double radiusHPlus[] = { 25, 45, 77 };

		double heightHMinus[] = { 24, 44, 58 };
		double radiusHMinus[] = { 29, 49, 63 };

		double height[] = heightHPlus;
		double radius[] = radiusHPlus;
		if (!isHadron) {
			height = heightElectron;
			radius = radiusElectron;
		}
		if (isHadron && !positive) {
			height = heightHMinus;
			radius = radiusHMinus;
		}
if(this.evtNumber==debugEvent && isHadron==false)
{
	System.out.println("check my event");
}
		double x[] = { part_DC_c1x[j], part_DC_c2x[j], part_DC_c3x[j] };
		double y[] = { part_DC_c1y[j], part_DC_c2y[j], part_DC_c3y[j] };

		int sec[] = { this.part_DC_sector[j] - 1, this.part_DC_sector[j] - 1, this.part_DC_sector[j] - 1 };
		for (int i = 0; i < 3; i++) {
			dcAllPos[i].fill(x[i], y[i]);
//this is a backwards rotation since the signs are flipped, accounting for -sin(x)=sin(-x), cos(x)=cos(-x)
			double x1_rot = y[i] * Math.sin(sec[i] * 60.0 * Pival / 180) + x[i] * Math.cos(sec[i] * 60.0 * Pival / 180);
			double y1_rot = y[i] * Math.cos(sec[i] * 60.0 * Pival / 180) - x[i] * Math.sin(sec[i] * 60.0 * Pival / 180);
			if(this.evtNumber==debugEvent && positive==false && isHadron==true)
			{
				//System.out.println("check my event");
				//System.out.println("det " +i + " x1 rot: "+ x1_rot +" y1_rot: "+y1_rot);
			}
				
			double slope = 1 / Math.tan(0.5 * angle * Pival / 180);
			double left = (height[i] - slope * y1_rot);
			double right = (height[i] + slope * y1_rot);
			double radius2 = Math.pow(radius[i], 2) - Math.pow(y1_rot, 2);
			if (x1_rot > left && x1_rot > right && Math.pow(x1_rot, 2) > radius2) {
				cut[i] = true;
				dcAcceptedPos[i].fill(x[i], y[i]);
			} else {
				cut[i] = false;	
			}
		}

		return (cut[0] && cut[1] && cut[2]);
	}

	boolean loadFTOF() {
		String bankName = "REC::Scintillator";
		if (!(m_event.hasBank(bankName))) {
			// System.out.println("couldn't find bank" + bankName);
			return false;
		} else {
			// System.out.println("got scintillator");
			HipoDataBank eventBank = (HipoDataBank) m_event.getBank(bankName);

			for (int current_part = 0; current_part < eventBank.rows(); current_part++) {
				// System.out.println("get pid");
				int pindex = eventBank.getInt("pindex", current_part);
				int sector = eventBank.getInt("sector", current_part);
				float time = eventBank.getFloat("time", current_part);
				int detID = eventBank.getInt("detector", current_part);
				float path = eventBank.getFloat("path", current_part);
				if (pindex >= maxArrSize) {
					System.out.println("too many particles is trajectory bank: " + pindex);
					break;
				}
				// entrance point to region 1, there are still multiple possible layers
				// however, as long as we associate the correct time and path, that should be
				// fine
				// also, we can just take the beta from the bank...
				if (detID == 12) {
					FTOFHit[pindex] = 1;
					FTOFSector[pindex] = sector;
					// System.out.println("setting ftof sector for " +pindex + " to "+ sector);
					FTOFTime[pindex] = time;
					FTOFPath[pindex] = path;
				} else {
					// this doesn't make sense
					// FTOFHit[pindex]=0;
					// FTOFSector[pindex]=-1;
					// FTOFTime[pindex]=-1;
				}

			}
		}
		return true;

	}

	int getHelicityAndSetRunEvtNumbers() {
		int helic=0;
		torus=(float)0.0;
		if(m_event.hasBank("RUN::config"))	
		{
			HipoDataBank runConfig = (HipoDataBank) m_event.getBank("RUN::config");
			torus=runConfig.getFloat("torus",0);
			solenoid=runConfig.getFloat("solenoid",0);
		}
		
		
		String bankName=new String();
		if(this.m_isMC)
			bankName="MC::Header";
		else
			bankName = "REC::Event";
		if (!m_event.hasBank(bankName)) {
			return 7;
		} 
		else {
			HipoDataBank eventBank = (HipoDataBank) m_event.getBank(bankName);
			if(this.m_isMC)
			{
			//	System.out.println("getting mc helicity");
				runNumber=eventBank.getInt("run",0);
				evtNumber=eventBank.getInt("event",0);
				float helicF= eventBank.getFloat("helicity", 0);	
				if(helicF<0)
					helic=-1;
				else
					helic=1;
				
			//System.out.println("helic is  " + helic);
			}
			else
			{
				runNumber=eventBank.getInt("NRUN",0);
				evtNumber=eventBank.getInt("NEVENT",0);
				helic= eventBank.getByte("Helic", 0);
				//System.out.println("helic is  " + helic);
			}
		}
		return helic;
	}

	boolean loadDCInfo() {
		String bankName = "REC::Traj";
		// needed to get DC sector
		String bankName2 = "TimeBasedTrkg::TBTracks";
		int idToSector[]=new int[this.maxArrSize];
		
		if (!useTimeBasedTracks)
			bankName2 = "REC::Track";
		if (!(m_event.hasBank(bankName) && m_event.hasBank(bankName2))) {
			// System.out.println("couldn't find bank" + bankName);
			return false;
		} else {
			// System.out.println("got trajectory or tracks");
			HipoDataBank eventBank = (HipoDataBank) m_event.getBank(bankName);
			HipoDataBank trkBank = (HipoDataBank) m_event.getBank(bankName2);

			for (int current_part = 0; current_part < trkBank.rows(); current_part++) {
				int sector = trkBank.getInt("sector", current_part);
				int trkPIndex=0;
				if(!useTimeBasedTracks)
				{
					 
					trkPIndex = trkBank.getInt("pindex", current_part);
					this.trkSectors[trkPIndex] = sector;
				}
				else
				{
					//time based tracks have only id, which we have to check against the 
					//id of the REC::Traj
					trkPIndex = trkBank.getInt("id", current_part);
					idToSector[trkPIndex]=sector;
				}
				
				// System.out.println("Trk " + current_part + " sector: "+sector);
			}

			for (int current_part = 0; current_part < eventBank.rows(); current_part++) {
				// System.out.println("get pid");
				// index in the track bank
				int index = eventBank.getInt("index", current_part);
				
				int pindex = eventBank.getInt("pindex", current_part);
				if(useTimeBasedTracks)
				{
					trkSectors[pindex]=idToSector[index];
				}
				int detID = eventBank.getInt("detId", current_part);

				if (pindex >= maxArrSize) {
					System.out.println("too many particles is trajectory bank: " + pindex);
					break;
				}

				part_DC_sector[pindex] = trkSectors[pindex];
				// System.out.println("got index " + index + " sector: " + trkSectors[index]);
				// entrance point to region 1
				if (detID == 12) {
					// trajectory doesn't have sector info...
					// int sector = eventBank.getInt("sector", current_part);
					// part_DC_sector[pindex]=sector;
					float x = eventBank.getFloat("x", current_part);
					part_DC_c1x[pindex] = x;
					float y = eventBank.getFloat("y", current_part);
					part_DC_c1y[pindex] = y;
				}
				if (detID == 24) {
					float x = eventBank.getFloat("x", current_part);
					part_DC_c2x[pindex] = x;
					float y = eventBank.getFloat("y", current_part);
					part_DC_c2y[pindex] = y;
				}
				if (detID == 36) {
					float x = eventBank.getFloat("x", current_part);
					part_DC_c3x[pindex] = x;
					float y = eventBank.getFloat("y", current_part);
					part_DC_c3y[pindex] = y;

					// counter++;
				}
			}
		}
		return true;

	}

	boolean loadCalInfo() {
		String bankName = "REC::Calorimeter";
		if (!(m_event.hasBank(bankName))) {
			// System.out.println("couldn't find bank" + bankName);
			return false;
		} else {
			// System.out.println("got calo");
			HipoDataBank eventBank = (HipoDataBank) m_event.getBank(bankName);
			int counter = 0;
			for (int current_part = 0; current_part < eventBank.rows(); current_part++) {
				// System.out.println("get pid");
				int pindex = eventBank.getInt("pindex", current_part);
				if (pindex >= maxArrSize) {
					System.out.println("too many particles is calorimeter bank: " + pindex);
					break;
				}
				int layer = eventBank.getInt("layer", current_part);
				// pcal is layer 1, 4 ECinner, 7 ECOuter
				float energy = eventBank.getFloat("energy", current_part);
				if (layer == 1) {
					// calPindex[counter]=pindex;
					int sector = eventBank.getInt("sector", current_part);
					part_Cal_PCAL_sector[pindex] = sector;
					float x = eventBank.getFloat("x", current_part);
					part_Cal_PCAL_x[pindex] = x;
					float y = eventBank.getFloat("y", current_part);
					part_Cal_PCAL_y[pindex] = y;
					
					float lu = eventBank.getFloat("lu", current_part);
					PCAL_lu[pindex] = lu;
					float lv = eventBank.getFloat("lv", current_part);
					PCAL_lv[pindex] = lv;
					float lw = eventBank.getFloat("lw", current_part);
					PCAL_lw[pindex] = lw;
					
					
					part_Cal_PCAL_E[pindex] = energy;
				}
				if (layer == 1 || layer == 4 || layer == 7) {
					// this was set to zero at the start of the event
					part_Cal_CalTotal_E[pindex] += energy;
				}
				counter++;
			}

		}
		return true;
	}

	boolean hasPCALInfo(int particleIndex) {
		String bankName = "REC::Calorimeter";

		if (!(m_event.hasBank(bankName))) {

			// System.out.println("couldn't find bank" + bankName);
			return false;
		} else {
			HipoDataBank eventBank = (HipoDataBank) m_event.getBank(bankName);
			for (int current_part = 0; current_part < eventBank.rows(); current_part++) {
				// System.out.println("get pid");
				int pindex = eventBank.getInt("pindex", current_part);
				if (pindex == particleIndex)
					return true;
			}
		}
		return true;
	}

	///////////
	/////// Stuff from Stefan for the

	// c) beta cuts

	boolean prot_beta_cut(int j, int run) {

		double prot_mean_p0[] = {1.00401, 1.00528, 1.00654, 1.00591,
				1.00188, 1.00275};
		double prot_mean_p1[] = {0.857477, 0.844641, 0.849807, 0.856546,
				0.868042, 0.865406};	
		double prot_sigma_p0[] = {0.000523482, 0.000246276, 0.000353689,
				0.0010485, 2.31291e-05, 0.000327357};
		double prot_sigma_p1[] = {0.0116411, 0.0127507, 0.0134067,
				0.0111433, 0.0115489, 0.0118189};

		double sigma_range = 3;

		double mean = 0;
		double sigma = 0;
		double upper_lim = 0;
		double lower_lim = 0;

		for (int k = 0; k < 6; k++) {
			if (this.FTOFSector[j] - 1 == k) {
				mean = prot_mean_p0[k] * part_p[j] / Math.sqrt(Math.pow(part_p[j], 2) + prot_mean_p1[k]);
				sigma = prot_sigma_p0[k] + prot_sigma_p1[k] / Math.sqrt(part_p[j]);
				upper_lim = mean + sigma_range * sigma;
				lower_lim = mean - sigma_range * sigma;
			}
		}

		if (this.part_beta[j] <= upper_lim && part_beta[j] >= lower_lim)
			return true;
		else
			return false;
	}

	boolean pip_beta_cut(int j, int run) {

		double pip_mean_p0[] = {0.999323, 0.999618, 0.999888, 0.999911,
				0.998791, 0.998889};
		double pip_mean_p1[] =  {0.0180815, 0.0170855, 0.017013, 0.0180826,
				0.0170883, 0.0163666};
		double pip_sigma_p0[] = {0.0010424, 0.00154015, 0.00179487,
				0.00224062, 0.000698921, 0.00103443};
		double pip_sigma_p1[] = {0.00666186, 0.00601759, 0.00596096,
				0.00509339, 0.00670074, 0.00668428};

		double sigma_range = 3;

		double mean = 0;
		double sigma = 0;
		double upper_lim = 0;
		double lower_lim = 0;

		for (int k = 0; k < 6; k++) {
			if (this.FTOFSector[j] - 1 == k) {
				mean = pip_mean_p0[k] * part_p[j] / Math.sqrt(Math.pow(part_p[j], 2) + pip_mean_p1[k]);
				sigma = pip_sigma_p0[k] + pip_sigma_p1[k] / Math.sqrt(part_p[j]);
				upper_lim = mean + sigma_range * sigma;
				lower_lim = mean - sigma_range * sigma;
			}
		}

		if (part_beta[j] <= upper_lim && part_beta[j] >= lower_lim)
			return true;
		else
			return false;
	}

	boolean pim_beta_cut(int j, int run) {

		double pim_mean_p0[] = {0.99633, 0.99649, 0.996569, 0.996175,
				0.99659, 0.996652};
		double pim_mean_p1[] =  {0.0152746, 0.0157846, 0.0153856, 0.0151136,
				0.014046, 0.0140012};
		double pim_sigma_p0[] = {6.31538e-05, -0.000321628, 0.000343613,
				3.57703e-05, -0.000221968, -0.000303001};
		double pim_sigma_p1[] = {0.00683594, 0.00730936, 0.00658419,
				0.00663091, 0.00704007, 0.0074271};

		double sigma_range = 3;

		double mean = 0;
		double sigma = 0;
		double upper_lim = 0;
		double lower_lim = 0;

		for (int k = 0; k < 6; k++) {
			if (FTOFSector[j] - 1 == k) {
				mean = pim_mean_p0[k] * part_p[j] / Math.sqrt(Math.pow(part_p[j], 2) + pim_mean_p1[k]);
				sigma = pim_sigma_p0[k] + pim_sigma_p1[k] / Math.sqrt(part_p[j]);
				upper_lim = mean + sigma_range * sigma;
				lower_lim = mean - sigma_range * sigma;
			}
		}

		if (part_beta[j] <= upper_lim && part_beta[j] >= lower_lim)
			return true;
		else
			return false;
	}

	boolean Kp_beta_cut(int j, int run) {

		double Kp_mean_p0[] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
		double Kp_mean_p1[] = { 0.24372, 0.24372, 0.24372, 0.24372, 0.24372, 0.24372 };
		double Kp_sigma_p0[] = { 0.00704184, 0.00112775, -8.45517e-05, 0.00346458, 0.00266644, 0.00314696 };
		double Kp_sigma_p1[] = { 0.00404813, 0.00984778, 0.0134949, 0.00843767, 0.00915269, 0.00921848 };

		double sigma_range = 3;

		double mean = 0;
		double sigma = 0;
		double upper_lim = 0;
		double lower_lim = 0;

		for (int k = 0; k < 6; k++) {
			if (FTOFSector[j] - 1 == k) {
				mean = Kp_mean_p0[k] * part_p[j] / Math.sqrt(Math.pow(part_p[j], 2) + Kp_mean_p1[k]);
				sigma = Kp_sigma_p0[k] + Kp_sigma_p1[k] / Math.sqrt(part_p[j]);
				upper_lim = mean + sigma_range * sigma;
				lower_lim = mean - sigma_range * sigma;
			}
		}

		if (part_beta[j] <= upper_lim && part_beta[j] >= lower_lim)
			return true;
		else
			return false;
	}

	boolean Km_beta_cut(int j, int run) {

		double Km_mean_p0[] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
		double Km_mean_p1[] = {0.244, 0.244, 0.244, 0.244, 0.244, 0.244};
		double Km_sigma_p0[] =  {0.0010424, 0.00154015, 0.00179487,
				0.00224062, 0.000698921, 0.00103443};
		double Km_sigma_p1[] = {0.00666186, 0.00601759, 0.00596096,
				0.00509339, 0.00670074, 0.00668428}; 

		double sigma_range = 3;

		double mean = 0;
		double sigma = 0;
		double upper_lim = 0;
		double lower_lim = 0;

		for (int k = 0; k < 6; k++) {
			if (FTOFSector[j] - 1 == k) {
				mean = Km_mean_p0[k] * part_p[j] / Math.sqrt(Math.pow(part_p[j], 2) + Km_mean_p1[k]);
				sigma = Km_sigma_p0[k] + Km_sigma_p1[k] / Math.sqrt(part_p[j]);
				upper_lim = mean + sigma_range * sigma;
				lower_lim = mean - sigma_range * sigma;
			}
		}

		if (part_beta[j] <= upper_lim && part_beta[j] >= lower_lim)
			return true;
		else
			return false;
	}

	boolean maximum_probability_cut(int j, int hypothesis, double conflvl, double anticonflvl) {

		// possible hypotheses which will be tested
		//
		// proton: 2212 pip: 211 Kp: 321
		// pim: -211 Km: -321
		//
		// particle variables:

		double sector = FTOFSector[j];
		double charge = part_charge[j];
		double mom = part_p[j];
		double beta = part_beta[j];

		// //////////////////////////////////////////////////////////////////////////////////////////////////
		// mean value and resolution for beta as a function of p for the different
		// sectors

		double prot_mean_p0[] = {1.00401, 1.00528, 1.00654, 1.00591,
				1.00188, 1.00275};
		 double prot_mean_p1[] = {0.857477, 0.844641, 0.849807, 0.856546,
				 0.868042, 0.865406};
				   double prot_sigma_p0[] = {0.000523482, 0.000246276, 0.000353689,
				 0.0010485, 2.31291e-05, 0.000327357};
				   double prot_sigma_p1[] = {0.0116411, 0.0127507, 0.0134067,
				 0.0111433, 0.0115489, 0.0118189};

				   double pip_mean_p0[] = {0.999323, 0.999618, 0.999888, 0.999911,
				 0.998791, 0.998889};
				   double pip_mean_p1[] = {0.0180815, 0.0170855, 0.017013, 0.0180826,
				 0.0170883, 0.0163666};
				   double pip_sigma_p0[] = {0.0010424, 0.00154015, 0.00179487,
				 0.00224062, 0.000698921, 0.00103443};
				   double pip_sigma_p1[] = {0.00666186, 0.00601759, 0.00596096,
				 0.00509339, 0.00670074, 0.00668428};

				   double Kp_mean_p0[] = {1.00, 1.00, 1.00, 1.00, 1.00, 1.00};
				 // literature
				   double Kp_mean_p1[] = {0.244, 0.244, 0.244, 0.244, 0.244, 0.244};
				 // literature
				   double Kp_sigma_p0[] = {0.0010424, 0.00154015, 0.00179487,
				 0.00224062, 0.000698921, 0.00103443};     // copied from pip
				   double Kp_sigma_p1[] = {0.00666186, 0.00601759, 0.00596096,
				 0.00509339, 0.00670074, 0.00668428};     // copied from pip

				   double pim_mean_p0[] = {0.99633, 0.99649, 0.996569, 0.996175,
				 0.99659, 0.996652};
				   double pim_mean_p1[] = {0.0152746, 0.0157846, 0.0153856, 0.0151136,
				 0.014046, 0.0140012};
				   double pim_sigma_p0[] = {6.31538e-05, -0.000321628, 0.000343613,
				 3.57703e-05, -0.000221968, -0.000303001};
				   double pim_sigma_p1[] = {0.00683594, 0.00730936, 0.00658419,
				 0.00663091, 0.00704007, 0.0074271};

				   double Km_mean_p0[] = {1.00, 1.00, 1.00, 1.00, 1.00, 1.00};
				 // literature
				   double Km_mean_p1[] = {0.244, 0.244, 0.244, 0.244, 0.244, 0.244};
				 // literature
				   double Km_sigma_p0[] = {6.31538e-05, -0.000321628, 0.000343613,
				 3.57703e-05, -0.000221968, -0.000303001};    // copied from pim
				   double Km_sigma_p1[] = {0.00683594, 0.00730936, 0.00658419,
				 0.00663091, 0.00704007, 0.0074271};              // copied from pim

		// //////////////////////////////////////////////////////////////////////////////////////////////////
		// population factors for the different particles (integrated over p)

		// intially no population weighting:

		double popfrac_proton = 1.0;
		double popfrac_pip = 1.0;
		double popfrac_Kp = 1.0;
		double popfrac_pim = 1.0;
		double popfrac_Km = 1.0;

		// momentum dependent population factor:

		if (population_weighting == true) {
			if (outbending == false) { // inbending (torus -1)
				popfrac_proton = -0.46270 + 1.497000000 * Math.pow(part_p[j], 1) - 0.9120000000 * Math.pow(part_p[j], 2)
						+ 0.242400000 * Math.pow(part_p[j], 3) - 0.0260900000 * Math.pow(part_p[j], 4)
						- 0.0009616000000 * Math.pow(part_p[j], 5) + 0.000520700 * Math.pow(part_p[j], 6)
						- 0.0000484100 * Math.pow(part_p[j], 7) + 0.000001559 * Math.pow(part_p[j], 8)
						+ 0.0000000071 * Math.pow(part_p[j], 9) - 0.0000000008185 * Math.pow(part_p[j], 10);
				popfrac_pip = 1.08600 - 0.739500000 * Math.pow(part_p[j], 1) + 0.3374000000 * Math.pow(part_p[j], 2)
						- 0.065160000 * Math.pow(part_p[j], 3) + 0.0060100000 * Math.pow(part_p[j], 4)
						- 0.0002440000000 * Math.pow(part_p[j], 5) + 0.000002776 * Math.pow(part_p[j], 6);
				popfrac_Kp = 0.08236 - 0.091010000 * Math.pow(part_p[j], 1) + 0.0784700000 * Math.pow(part_p[j], 2)
						- 0.021870000 * Math.pow(part_p[j], 3) + 0.0023840000 * Math.pow(part_p[j], 4)
						- 0.0000307300000 * Math.pow(part_p[j], 5) - 0.000011330 * Math.pow(part_p[j], 6)
						+ 0.0000005536 * Math.pow(part_p[j], 7);

				popfrac_pim = 0.95340 + 0.021960000 * Math.pow(part_p[j], 1) - 0.0405300000 * Math.pow(part_p[j], 2)
						+ 0.009071000 * Math.pow(part_p[j], 3) - 0.0006846000 * Math.pow(part_p[j], 4)
						+ 0.0000155900000 * Math.pow(part_p[j], 5);
				popfrac_Km = -0.18360 + 0.714600000 * Math.pow(part_p[j], 1) - 0.7519000000 * Math.pow(part_p[j], 2)
						+ 0.390900000 * Math.pow(part_p[j], 3) - 0.1085000000 * Math.pow(part_p[j], 4)
						+ 0.0170200000000 * Math.pow(part_p[j], 5) - 0.001513000 * Math.pow(part_p[j], 6)
						+ 0.0000708600 * Math.pow(part_p[j], 7) - 0.000001348 * Math.pow(part_p[j], 8);
			}
			if (outbending == true) { // outbending (torus +1)
				popfrac_proton = -0.43870 + 1.5410000000 * Math.pow(part_p[j], 1)
						- 0.21180000000 * Math.pow(part_p[j], 2) - 1.0700000 * Math.pow(part_p[j], 3)
						+ 0.95490000 * Math.pow(part_p[j], 4) - 0.38710000 * Math.pow(part_p[j], 5)
						+ 0.0858300000 * Math.pow(part_p[j], 6) - 0.00966400000 * Math.pow(part_p[j], 7)
						+ 0.0001328 * Math.pow(part_p[j], 8) + 0.0001076 * Math.pow(part_p[j], 9)
						- 0.00001395 * Math.pow(part_p[j], 10	) + 0.0000007545 * Math.pow(part_p[j], 11)
						- 0.00000001582 * Math.pow(part_p[j], 12);
				popfrac_pip = 0.90020 - 0.7225000000 * Math.pow(part_p[j], 1) + 0.41450000000 * Math.pow(part_p[j], 2)
						- 0.1149000 * Math.pow(part_p[j], 3) + 0.01728000 * Math.pow(part_p[j], 4)
						- 0.00132500 * Math.pow(part_p[j], 5) + 0.0000400800 * Math.pow(part_p[j], 6);
				popfrac_Kp = 0.17180 - 0.1994000000 * Math.pow(part_p[j], 1) + 0.11450000000 * Math.pow(part_p[j], 2)
						- 0.0190200 * Math.pow(part_p[j], 3) + 0.00040970 * Math.pow(part_p[j], 4)
						+ 0.00012380 * Math.pow(part_p[j], 5) - 0.0000070660 * Math.pow(part_p[j], 6);

				popfrac_pim = 0.95950 + 0.0063400000 * Math.pow(part_p[j], 1) - 0.030890000000 * Math.pow(part_p[j], 2)
						+ 0.0052940 * Math.pow(part_p[j], 3) - 0.00024860 * Math.pow(part_p[j], 4);
				popfrac_Km = 0.02474 + 0.0132600000 * Math.pow(part_p[j], 1) + 0.024240000000 * Math.pow(part_p[j], 2)
						- 0.0045600 * Math.pow(part_p[j], 3) + 0.00022150 * Math.pow(part_p[j], 4);
			}
		}

		if (charge > 0) {
			if (hypothesis <= 0)
				return false; // charge does not match with hypothesis

			for (int k = 0; k < 6; k++) {
				if (sector - 1 == k) {

					double mean_prot = prot_mean_p0[k] * part_p[j]
							/ Math.sqrt(Math.pow(part_p[j], 2) + prot_mean_p1[k]);
					double sigma_prot = prot_sigma_p0[k] + prot_sigma_p1[k] / Math.sqrt(part_p[j]);
					double prob_prot = popfrac_proton * (1 / (sigma_prot * Math.sqrt(2 * 3.14159)))
							* Math.exp(-0.5 * Math.pow((beta - mean_prot) / sigma_prot, 2));
					double conf_prot = 100 * (1.0 - org.apache.commons.math3.special.Erf
							.erf(Math.abs(beta - mean_prot) / sigma_prot / Math.sqrt(2.0)));

					double mean_pip = pip_mean_p0[k] * part_p[j] / Math.sqrt(Math.pow(part_p[j], 2) + pip_mean_p1[k]);
					double sigma_pip = pip_sigma_p0[k] + pip_sigma_p1[k] / Math.sqrt(part_p[j]);
					double prob_pip = popfrac_pip * (1 / (sigma_pip * Math.sqrt(2 * 3.14159)))
							* Math.exp(-0.5 * Math.pow((beta - mean_pip) / sigma_pip, 2));
					double conf_pip = 100 * (1.0 - org.apache.commons.math3.special.Erf
							.erf(Math.abs(beta - mean_pip) / sigma_pip / Math.sqrt(2.0)));

					// double mean_Kp = Kp_mean_p0[k] * part_p[j] / Math.sqrt(Math.pow(part_p[j],2)
					// + Kp_mean_p1[k]);
					double mean_Kp = part_p[j] / Math.sqrt(Math.pow(part_p[j], 2) + Math.pow(0.493677, 2));
					double sigma_Kp = Kp_sigma_p0[k] + Kp_sigma_p1[k] / Math.sqrt(part_p[j]);
					double prob_Kp = popfrac_Kp * (1 / (sigma_Kp * Math.sqrt(2 * 3.14159)))
							* Math.exp(-0.5 * Math.pow((beta - mean_Kp) / sigma_Kp, 2));
					double conf_Kp = 100 * (1.0 - org.apache.commons.math3.special.Erf
							.erf(Math.abs(beta - mean_Kp) / sigma_Kp / Math.sqrt(2.0)));
					// prob_Kp = 0; // overwrite Kaons
					// conf_Kp = 0; // overwrite Kaons

					if (prob_prot > prob_pip && prob_prot > prob_Kp && hypothesis == 2212 && conf_prot > conflvl
							&& conf_pip < anticonflvl && conf_Kp < anticonflvl) {
						// cout << "proton - probability: " << prob_prot << " confidence: " << conf_prot
						// << " ( mom: " << mom << " beta:" << beta << " mean: " << mean_prot << "
						// sigma: " << sigma_prot << " )" << endl;
						// cout << "pip - probability: " << prob_pip << " confidence: " << conf_pip << "
						// ( mom: " << mom << " beta:" << beta << " mean: " << mean_pip << " sigma: " <<
						// sigma_pip << " )" << endl;
						// cout << "Kp - probability: " << prob_Kp << " confidence: " << conf_Kp << " (
						// mom: " << mom << " beta:" << beta << " mean: " << mean_Kp << " sigma: " <<
						// sigma_Kp << " )" << endl;
						// cout << endl;
						return true;
					}
					if (prob_pip > prob_prot && prob_pip > prob_Kp && hypothesis == 211 && conf_pip > conflvl
							&& conf_prot < anticonflvl && conf_Kp < anticonflvl) {
						// cout << "proton - probability: " << prob_prot << " confidence: " << conf_prot
						// << " ( mom: " << mom << " beta:" << beta << " mean: " << mean_prot << "
						// sigma: " << sigma_prot << " )" << endl;
						// cout << "pip - probability: " << prob_pip << " confidence: " << conf_pip << "
						// ( mom: " << mom << " beta:" << beta << " mean: " << mean_pip << " sigma: " <<
						// sigma_pip << " )" << endl;
						// cout << "Kp - probability: " << prob_Kp << " confidence: " << conf_Kp << " (
						// mom: " << mom << " beta:" << beta << " mean: " << mean_Kp << " sigma: " <<
						// sigma_Kp << " )" << endl;
						// cout << endl;
						return true;
					}
					if (prob_Kp > prob_prot && prob_Kp > prob_pip && hypothesis == 321 && conf_Kp > conflvl
							&& conf_prot < anticonflvl && conf_pip < anticonflvl) {
						// cout << "proton - probability: " << prob_prot << " confidence: " << conf_prot
						// << " ( mom: " << mom << " beta:" << beta << " mean: " << mean_prot << "
						// sigma: " << sigma_prot << " )" << endl;
						// cout << "pip - probability: " << prob_pip << " confidence: " << conf_pip << "
						// ( mom: " << mom << " beta:" << beta << " mean: " << mean_pip << " sigma: " <<
						// sigma_pip << " )" << endl;
						// cout << "Kp - probability: " << prob_Kp << " confidence: " << conf_Kp << " (
						// mom: " << mom << " beta:" << beta << " mean: " << mean_Kp << " sigma: " <<
						// sigma_Kp << " )" << endl;
						// cout << endl;
						return true;
					}
				}
			}
		}

		if (charge < 0) {
			if (hypothesis >= 0)
				return false; // charge does not match with hypothesis

			for (int k = 0; k < 6; k++) {
				if (sector - 1 == k) {

					double mean_pim = pim_mean_p0[k] * part_p[j] / Math.sqrt(Math.pow(part_p[j], 2) + pim_mean_p1[k]);
					double sigma_pim = pim_sigma_p0[k] + pim_sigma_p1[k] / Math.sqrt(part_p[j]);
					double prob_pim = popfrac_pim * (1 / (sigma_pim * Math.sqrt(2 * 3.14159)))
							* Math.exp(-0.5 * Math.pow((beta - mean_pim) / sigma_pim, 2));
					double conf_pim = 100 * (1.0 - org.apache.commons.math3.special.Erf
							.erf(Math.abs(beta - mean_pim) / sigma_pim / Math.sqrt(2.0)));

					// double mean_Km = Km_mean_p0[k] * part_p[j] / Math.sqrt(Math.pow(part_p[j],2)
					// + Km_mean_p1[k]);
					double mean_Km = part_p[j] / Math.sqrt(Math.pow(part_p[j], 2) + Math.pow(0.493677, 2));
					double sigma_Km = Km_sigma_p0[k] + Km_sigma_p1[k] / Math.sqrt(part_p[j]);
					double prob_Km = popfrac_Km * (1 / (sigma_Km * Math.sqrt(2 * 3.14159)))
							* Math.exp(-0.5 * Math.pow((beta - mean_Km) / sigma_Km, 2));
					double conf_Km = 100 * (1.0 - org.apache.commons.math3.special.Erf
							.erf(Math.abs(beta - mean_Km) / sigma_Km / Math.sqrt(2.0)));
					// prob_Km = 0; // overwrite Kaons
					// conf_Km = 0; // overwrite Kaons

					if (prob_pim > prob_Km && hypothesis == -211 && conf_pim > conflvl && conf_Km < anticonflvl) {
						// cout << "pim - probability: " << prob_pim << " confidence: " << conf_pim << "
						// ( mom: " << mom << " beta:" << beta << " mean: " << mean_pim << " sigma: " <<
						// sigma_pim << " )" << endl;
						// cout << "Km - probability: " << prob_Km << " confidence: " << conf_Km << " (
						// mom: " << mom << " beta:" << beta << " mean: " << mean_Km << " sigma: " <<
						// sigma_Km << " )" << endl;
						// cout << endl;
						return true;
					}
					if (prob_Km > prob_pim && hypothesis == -321 && conf_Km > conflvl && conf_pim < anticonflvl) {
						// cout << "pim - probability: " << prob_pim << " confidence: " << conf_pim << "
						// ( mom: " << mom << " beta:" << beta << " mean: " << mean_pim << " sigma: " <<
						// sigma_pim << " )" << endl;
						// cout << "Km - probability: " << prob_Km << " confidence: " << conf_Km << " (
						// mom: " << mom << " beta:" << beta << " mean: " << mean_Km << " sigma: " <<
						// sigma_Km << " )" << endl;
						// cout << endl;
						return true;
					}
				}
			}
		}

		return false; // return false if no particle is clearly identified or if the charge is 0.
	}

	// g) delta vz cuts

	boolean prot_delta_vz_cut(int j) {

		double mean = 0.4543;
		double sigma = 4.558;
		double dvz_min = mean - 3 * sigma;
		double dvz_max = mean + 3 * sigma;

		if (part_vz[j] > dvz_min && part_vz[j] < dvz_max)
			return true;
		else
			return false;
	}

	boolean pip_delta_vz_cut(int j) {

		double mean = 0.5252;
		double sigma = 4.483;
		double dvz_min = mean - 3 * sigma;
		double dvz_max = mean + 3 * sigma;

		if (part_vz[j] > dvz_min && part_vz[j] < dvz_max)
			return true;
		else
			return false;
	}

	boolean pim_delta_vz_cut(int j) {

		double mean = 0.54;
		double sigma = 4.699;
		double dvz_min = mean - 3 * sigma;
		double dvz_max = mean + 3 * sigma;

		if (part_vz[j] > dvz_min && part_vz[j] < dvz_max)
			return true;
		else
			return false;
	}

	boolean Kp_delta_vz_cut(int j) {

		double mean = 0.5;
		double sigma = 4.429;
		double dvz_min = mean - 3 * sigma;
		double dvz_max = mean + 3 * sigma;

		if (part_vz[j] > dvz_min && part_vz[j] < dvz_max)
			return true;
		else
			return false;
	}

	boolean Km_delta_vz_cut(int j) {

		double mean = 0.06038;
		double sigma = 3.984;
		double dvz_min = mean - 3 * sigma;
		double dvz_max = mean + 3 * sigma;

		if (part_vz[j] > dvz_min && part_vz[j] < dvz_max)
			return true;
		else
			return false;
	}

	boolean CD_maximum_probability_cut(int j, int hypothesis, double conflvl, double anticonflvl, int run) {

		// possible hypotheses which will be tested
		//
		// proton: 2212 pip: 211 Kp: 321
		// pim: -211 Km: -321
		//
		// particle variables:

		double charge = part_charge[j];
		double mom = part_p[j];
		double beta = part_beta[j];

		// //////////////////////////////////////////////////////////////////////////////////////////////////
		// mean value and resolution for beta as a function of p for the different
		// sectors

		double prot_mean_p0 = -0.0157773;
		  double prot_mean_p1 = 0.992786;
		  double prot_mean_p2 = 0.886921;
		  double prot_sigma_p0 = 0.129247;
		  double prot_sigma_p1 = -0.0606478;

		  double pip_mean_p0 = -1.40262;
		  double pip_mean_p1 = 2.37572;
		  double pip_mean_p2 = 0.0100619;
		  double pip_sigma_p0 = 0.11545;
		  double pip_sigma_p1 = -0.00569233;

		  double pim_mean_p0 = -0.309342;
		  double pim_mean_p1 = 1.28698;
		  double pim_mean_p2 = 0.0228848;
		  double pim_sigma_p0 = 0.0880416;
		  double pim_sigma_p1 = 0.00147085;

		  double Kp_mean_p0 = 0.00000;
		  double Kp_mean_p1 = 1.00000;
		  double Kp_mean_p2 = 0.24372;         // lit. Kaon mass
		  double Kp_sigma_p0 = 0.11545;        // copied from pip
		  double Kp_sigma_p1 = -0.00569233;    // copied from pip


		  double Km_mean_p0 = 0.00000;
		  double Km_mean_p1 = 1.00000;
		  double Km_mean_p2 = 0.24372;        // lit. Kaon mass
		  double Km_sigma_p0 = 0.0880416;     // copied from pim
		  double Km_sigma_p1 = 0.00147085;    // copied from pim

		// //////////////////////////////////////////////////////////////////////////////////////////////////
		// population factors for the different particles (integrated over p)

		// intially no population weighting:

		double popfrac_proton = 1.0;
		double popfrac_pip = 1.0;
		double popfrac_Kp = 1.0;
		double popfrac_pim = 1.0;
		double popfrac_Km = 1.0;

		// momentum dependent population factor:

		if (population_weighting_CD == true) {
			if (outbending == false) { // inbending (torus -1)
				popfrac_proton = 1.27200 - 6.6060000000 * Math.pow(part_p[j], 1)
						+ 12.8100000000 * Math.pow(part_p[j], 2) - 10.61000 * Math.pow(part_p[j], 3)
						+ 4.503000000 * Math.pow(part_p[j], 4) - 1.034000000 * Math.pow(part_p[j], 5)
						+ 0.1225000000 * Math.pow(part_p[j], 6) - 0.00588400000 * Math.pow(part_p[j], 7);
				if (part_p[j] > 5.0)
					popfrac_proton = 0.001;
				popfrac_pip = 0.87680 + 1.1150000000 * Math.pow(part_p[j], 1) - 3.94000000000 * Math.pow(part_p[j], 2)
						+ 3.679000 * Math.pow(part_p[j], 3) - 1.659000000 * Math.pow(part_p[j], 4)
						+ 0.425800000 * Math.pow(part_p[j], 5) - 0.0654700000 * Math.pow(part_p[j], 6)
						+ 0.00597200000 * Math.pow(part_p[j], 7) - 0.000298 * Math.pow(part_p[j], 8)
						+ 0.000006271 * Math.pow(part_p[j], 9);
				if (part_p[j] > 6.0)
					popfrac_pip = 0.998;
				popfrac_Kp = 0.04908 + 0.0908700000 * Math.pow(part_p[j], 1) - 0.032710000000 * Math.pow(part_p[j], 2)
						+ 0.002760 * Math.pow(part_p[j], 3);
				if (part_p[j] > 5.5)
					popfrac_Kp = 0.001;
				popfrac_pim = 1.003000 - 0.324900000 * Math.pow(part_p[j], 1) + 0.184200000 * Math.pow(part_p[j], 2)
						- 0.0335700000 * Math.pow(part_p[j], 3) + 0.00164000 * Math.pow(part_p[j], 4)
						+ 0.0001136 * Math.pow(part_p[j], 5) - 0.000009368 * Math.pow(part_p[j], 6);
				if (part_p[j] > 3.5)
					popfrac_pim = 0.999;
				popfrac_Km = 0.065570 - 0.155300000 * Math.pow(part_p[j], 1) + 0.652600000 * Math.pow(part_p[j], 2)
						- 0.5292000000 * Math.pow(part_p[j], 3) + 0.15800000 * Math.pow(part_p[j], 4)
						- 0.0162400 * Math.pow(part_p[j], 5);
				if (part_p[j] > 3.5)
					popfrac_Km = 0.001;
			}
			if (outbending == true) { // outbending (torus +1)
				popfrac_proton = 0.091220 - 0.774600000 * Math.pow(part_p[j], 1) + 3.121000000 * Math.pow(part_p[j], 2)
						- 3.1810000000 * Math.pow(part_p[j], 3) + 1.59300000 * Math.pow(part_p[j], 4)
						- 0.4622000 * Math.pow(part_p[j], 5) + 0.081690000 * Math.pow(part_p[j], 6)
						- 0.008687000 * Math.pow(part_p[j], 7) + 0.0005144000 * Math.pow(part_p[j], 8)
						- 0.00001281 * Math.pow(part_p[j], 9);
				if (part_p[j] > 8.0)
					popfrac_proton = 0.001;
				popfrac_pip = 1.223000 - 1.400000000 * Math.pow(part_p[j], 1) + 0.708100000 * Math.pow(part_p[j], 2)
						- 0.1232000000 * Math.pow(part_p[j], 3) + 0.00514000 * Math.pow(part_p[j], 4)
						+ 0.0005800 * Math.pow(part_p[j], 5) - 0.000013630 * Math.pow(part_p[j], 6)
						- 0.000006621 * Math.pow(part_p[j], 7) + 0.0000003665 * Math.pow(part_p[j], 8);
				if (part_p[j] > 7.0)
					popfrac_pip = 0.998;
				popfrac_Kp = -0.06797 + 0.294800000 * Math.pow(part_p[j], 1) - 0.137800000 * Math.pow(part_p[j], 2)
						+ 0.0231100000 * Math.pow(part_p[j], 3) - 0.00133500 * Math.pow(part_p[j], 4);
				if (part_p[j] > 6.0)
					popfrac_Kp = 0.001;
				popfrac_pim = 0.93350 + 0.1171000000 * Math.pow(part_p[j], 1) - 0.23500000000 * Math.pow(part_p[j], 2)
						- 0.254100 * Math.pow(part_p[j], 3) + 0.464600000 * Math.pow(part_p[j], 4)
						- 0.255500000 * Math.pow(part_p[j], 5) + 0.0702800000 * Math.pow(part_p[j], 6)
						- 0.01024000000 * Math.pow(part_p[j], 7) + 0.000606 * Math.pow(part_p[j], 8)
						+ 0.000034870 * Math.pow(part_p[j], 9) - 0.000008172 * Math.pow(part_p[j], 10)
						+ 0.0000005102 * Math.pow(part_p[j], 11) - 0.00000001146 * Math.pow(part_p[j], 12);
				if (part_p[j] > 4.0)
					popfrac_pim = 0.999;
				popfrac_Km = 0.10210 - 0.2916000000 * Math.pow(part_p[j], 1) + 0.72860000000 * Math.pow(part_p[j], 2)
						- 0.486800 * Math.pow(part_p[j], 3) + 0.124600000 * Math.pow(part_p[j], 4)
						- 0.011080000 * Math.pow(part_p[j], 5);
				if (part_p[j] > 4.0)
					popfrac_Km = 0.001;
			}
		}

		if (charge > 0) {

			if (hypothesis <= 0)
				return false; // charge does not match with hypothesis

			double mean_prot = prot_mean_p0
					+ prot_mean_p1 * part_p[j] / Math.sqrt(Math.pow(part_p[j], 2) + prot_mean_p2);
			double sigma_prot = prot_sigma_p0 + prot_sigma_p1 / Math.sqrt(part_p[j]);
			double prob_prot = popfrac_proton * (1 / (sigma_prot * Math.sqrt(2 * 3.14159)))
					* Math.exp(-0.5 * Math.pow((beta - mean_prot) / sigma_prot, 2));
			double conf_prot = 100 * (1.0 - org.apache.commons.math3.special.Erf
					.erf(Math.abs(beta - mean_prot) / sigma_prot / Math.sqrt(2.0)));

			double mean_pip = pip_mean_p0 + pip_mean_p1 * part_p[j] / Math.sqrt(Math.pow(part_p[j], 2) + pip_mean_p2);
			double sigma_pip = pip_sigma_p0 + pip_sigma_p1 / Math.sqrt(part_p[j]);
			double prob_pip = popfrac_pip * (1 / (sigma_pip * Math.sqrt(2 * 3.14159)))
					* Math.exp(-0.5 * Math.pow((beta - mean_pip) / sigma_pip, 2));
			double conf_pip = 100 * (1.0
					- org.apache.commons.math3.special.Erf.erf(Math.abs(beta - mean_pip) / sigma_pip / Math.sqrt(2.0)));

			double mean_Kp = Kp_mean_p0 + Kp_mean_p1 * part_p[j] / Math.sqrt(Math.pow(part_p[j], 2) + Kp_mean_p2);
			double sigma_Kp = Kp_sigma_p0 + Kp_sigma_p1 / Math.sqrt(part_p[j]);
			double prob_Kp = popfrac_Kp * (1 / (sigma_Kp * Math.sqrt(2 * 3.14159)))
					* Math.exp(-0.5 * Math.pow((beta - mean_Kp) / sigma_Kp, 2));
			double conf_Kp = 100 * (1.0
					- org.apache.commons.math3.special.Erf.erf(Math.abs(beta - mean_Kp) / sigma_Kp / Math.sqrt(2.0)));

			if (prob_prot > prob_pip && prob_prot > prob_Kp && hypothesis == 2212 && conf_prot > conflvl
					&& conf_pip < anticonflvl && conf_Kp < anticonflvl && part_p[j] > 0.3) {
				// cout << "proton - probability: " << prob_prot << " confidence: " << conf_prot
				// << " ( mom: " << mom << " beta:" << beta << " mean: " << mean_prot << "
				// sigma: " << sigma_prot << " )" << endl;
				// cout << "pip - probability: " << prob_pip << " confidence: " << conf_pip << "
				// ( mom: " << mom << " beta:" << beta << " mean: " << mean_pip << " sigma: " <<
				// sigma_pip << " )" << endl;
				// cout << "Kp - probability: " << prob_Kp << " confidence: " << conf_Kp << " (
				// mom: " << mom << " beta:" << beta << " mean: " << mean_Kp << " sigma: " <<
				// sigma_Kp << " )" << endl;
				// cout << endl;
				return true;
			}
			if (prob_pip > prob_prot && prob_pip > prob_Kp && hypothesis == 211 && conf_pip > conflvl
					&& conf_prot < anticonflvl && conf_Kp < anticonflvl && part_p[j] > 0.7) {
				// cout << "proton - probability: " << prob_prot << " confidence: " << conf_prot
				// << " ( mom: " << mom << " beta:" << beta << " mean: " << mean_prot << "
				// sigma: " << sigma_prot << " )" << endl;
				// cout << "pip - probability: " << prob_pip << " confidence: " << conf_pip << "
				// ( mom: " << mom << " beta:" << beta << " mean: " << mean_pip << " sigma: " <<
				// sigma_pip << " )" << endl;
				// cout << "Kp - probability: " << prob_Kp << " confidence: " << conf_Kp << " (
				// mom: " << mom << " beta:" << beta << " mean: " << mean_Kp << " sigma: " <<
				// sigma_Kp << " )" << endl;
				// cout << endl;
				return true;
			}
			if (prob_Kp > prob_prot && prob_Kp > prob_pip && hypothesis == 321 && conf_Kp > conflvl
					&& conf_prot < anticonflvl && conf_pip < anticonflvl && part_p[j] > 0.7) {
				// cout << "proton - probability: " << prob_prot << " confidence: " << conf_prot
				// << " ( mom: " << mom << " beta:" << beta << " mean: " << mean_prot << "
				// sigma: " << sigma_prot << " )" << endl;
				// cout << "pip - probability: " << prob_pip << " confidence: " << conf_pip << "
				// ( mom: " << mom << " beta:" << beta << " mean: " << mean_pip << " sigma: " <<
				// sigma_pip << " )" << endl;
				// cout << "Kp - probability: " << prob_Kp << " confidence: " << conf_Kp << " (
				// mom: " << mom << " beta:" << beta << " mean: " << mean_Kp << " sigma: " <<
				// sigma_Kp << " )" << endl;
				// cout << endl;
				return true;
			}
		}

		if (charge < 0) {

			if (hypothesis >= 0)
				return false; // charge does not match with hypothesis

			double mean_pim = pim_mean_p0 + pim_mean_p1 * part_p[j] / Math.sqrt(Math.pow(part_p[j], 2) + pim_mean_p2);
			double sigma_pim = pim_sigma_p0 + pim_sigma_p1 / Math.sqrt(part_p[j]);
			double prob_pim = popfrac_pim * (1 / (sigma_pim * Math.sqrt(2 * 3.14159)))
					* Math.exp(-0.5 * Math.pow((beta - mean_pim) / sigma_pim, 2));
			double conf_pim = 100 * (1.0 - (Math.abs(beta - mean_pim) / sigma_pim / Math.sqrt(2.0)));

			double mean_Km = Km_mean_p0 + Km_mean_p1 * part_p[j] / Math.sqrt(Math.pow(part_p[j], 2) + Km_mean_p2);
			double sigma_Km = Km_sigma_p0 + Km_sigma_p1 / Math.sqrt(part_p[j]);
			double prob_Km = popfrac_Km * (1 / (sigma_Km * Math.sqrt(2 * 3.14159)))
					* Math.exp(-0.5 * Math.pow((beta - mean_Km) / sigma_Km, 2));
			double conf_Km = 100 * (1.0
					- org.apache.commons.math3.special.Erf.erf(Math.abs(beta - mean_Km) / sigma_Km / Math.sqrt(2.0)));

			if (prob_pim > prob_Km && hypothesis == -211 && conf_pim > conflvl && conf_Km < anticonflvl
					&& part_p[j] > 0.7) {
				// cout << "pim - probability: " << prob_pim << " confidence: " << conf_pim << "
				// ( mom: " << mom << " beta:" << beta << " mean: " << mean_pim << " sigma: " <<
				// sigma_pim << " )" << endl;
				// cout << "Km - probability: " << prob_Km << " confidence: " << conf_Km << " (
				// mom: " << mom << " beta:" << beta << " mean: " << mean_Km << " sigma: " <<
				// sigma_Km << " )" << endl;
				// cout << endl;
				return true;
			}
			if (prob_Km > prob_pim && hypothesis == -321 && conf_Km > conflvl && conf_pim < anticonflvl
					&& part_p[j] > 0.7) {
				// cout << "pim - probability: " << prob_pim << " confidence: " << conf_pim << "
				// ( mom: " << mom << " beta:" << beta << " mean: " << mean_pim << " sigma: " <<
				// sigma_pim << " )" << endl;
				// cout << "Km - probability: " << prob_Km << " confidence: " << conf_Km << " (
				// mom: " << mom << " beta:" << beta << " mean: " << mean_Km << " sigma: " <<
				// sigma_Km << " )" << endl;
				// cout << endl;
				return true;
			}
		}

		return false; // return false if no particle is clearly identified or if the charge is 0.
	}

}