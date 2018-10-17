package NovelFitters;

import java.io.File;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.nio.file.PathMatcher;
import java.nio.file.Paths;

import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.io.hipo.HipoDataEvent;
import org.jlab.io.hipo.HipoDataSource;
import NovelFitters.NovelBaseFitter;
import NovelFitters.MyParticle;
import org.jlab.io.hipo.*;
import org.jlab.io.base.DataEvent;

import java.util.ArrayList;
import java.util.List;

import org.jlab.clas.physics.*;
import org.jlab.groot.data.H2F;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.groot.fitter.ParallelSliceFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.data.H1F;

import java.io.*;
import java.nio.file.*;
import java.nio.file.attribute.*;
import static java.nio.file.FileVisitResult.*;
import static java.nio.file.FileVisitOption.*;
import java.util.*;

public class FitterTester {
	protected NovelBaseFitter novel_fitter;
	protected NovelBaseFitter novel_fitterMC;
	public boolean isMC;
	public int numMatch;
	public int numAllPairs;
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		if (args.length == 0) {
			// exits program if input directory not specified
			System.out.println("ERROR: Please enter a hipo file as the first argument");
			System.exit(0);
		}
		FitterTester tester = new FitterTester();
		tester.isMC=false;
		//analyzer.isMC=false;
		tester.analyze(args);
		System.out.println("num match: " + tester.numMatch + " all pairs: "+ tester.numAllPairs + " fraction: " + (float)tester.numMatch/(float)tester.numAllPairs);
	}

	
	public void analyze(String[] args) {
		
		NovelBaseFitter.useStefanElectronCuts=true;
		NovelBaseFitter.useStefanPIDCuts=false;
		NovelBaseFitter.useStefanHadronCuts=true;
		NovelBaseFitter.useTimeBasedTracks=false;
		NovelBaseFitter.noVertexCut=true;
		HipoDataSource reader = new HipoDataSource();
		// define fitter class, argument is the beam energy (GeV)
		 novel_fitter = new NovelBaseFitter(10.6,false,false);
		 
		//novel_fitter = new NovelBaseFitter(10.6, true, true);
		 if(isMC)
			 novel_fitterMC = new NovelBaseFitter(10.6, true, true);
		// define filter to apply to your events
		// here we look for events with one electron (11), one photon (22) (change to no
		// photon) and any number of other
		// positively charged particles (X+), negatively charged particles (X-) or
		// neutral
		// particles (Xn)
		//EventFilter filter = new EventFilter("Xn");
		//EventFilter filter = new EventFilter("11:X+:X-:Xn");
		EventFilter filter = new EventFilter("11:+211:-211:X+:X-:Xn");
		File folder = new File(args[0]);
		File[] listOfFiles = folder.listFiles();
		for (int iF = 0; iF < listOfFiles.length; iF++) {
			if (listOfFiles[iF].isFile()) {
				System.out.println("File " + listOfFiles[iF].getName());
				PathMatcher matcher = FileSystems.getDefault().getPathMatcher("glob:*.{hipo}");

				Path filename = Paths.get(listOfFiles[iF].getName());
				if (matcher.matches(filename)) {
					//one asymmetry file per data input file
					
					System.out.println("matched" + filename);

					reader.open(args[0] + listOfFiles[iF].getName()); // open hipo file

					while (reader.hasEvent() == true) { // cycle through events
						// load next event in the hipo file
					//	System.out.println("new event----\n\n");
						HipoDataEvent event = (HipoDataEvent) reader.getNextEvent();
						
						// apply fitter to your data event
						PhysicsEvent generic_Event = novel_fitter.getPhysicsEvent(event);
						PhysicsEvent generic_EventMC=new PhysicsEvent();
						if(isMC)
						{
							
								generic_EventMC = novel_fitterMC.getPhysicsEvent(event);
						
						if(generic_Event.count()>1)
						{
					//	System.out.println("we have " + generic_Event.count() + " rec events and "+ generic_EventMC.count() + " mc events");
						}
							//-->
							associateMCWithData(generic_Event, generic_EventMC);
						
						}
						doDiHadrons(generic_Event,generic_EventMC,novel_fitter,novel_fitterMC);
						
						
					//	System.out.println("helicity is "+ generic_Event.)
						// novel_fitter.Walt);
						//51831184
						if(generic_Event.count()>1)
						{
						//	System.out.println("looking at event with " + generic_Event.count() + " particles ");
						}
						if (filter.isValid(generic_Event) == true) { // apply filter to current event
							// look at all particles
							if(novel_fitter.getEvtNumber()==novel_fitter.debugEvent)
							{
								System.out.println("event " + novel_fitter.debugEvent+" filtered event");
								System.out.println("we have " +generic_Event.count());
							}
							
							for (int i = 0; i < generic_Event.count(); i++) 	
							{
								MyParticle part = (MyParticle) generic_Event.getParticle(i);
								int sec= part.FTOFsector;
								if(sec<=0)
									sec=6;
								
								//System.out.println("sec: " + sec+ ", beta is " + part.beta);
								if(part.beta<10.0 && part.beta>-10.0)
								{
									
									if(sec!=6)
									{
										//System.out.println("fill with beta "+ part.beta);
									}
								}
								//System.out.println("time is: " + part.FTOFTime + " sector: " + part.FTOFsector);
								// System.out.println("matching mc particle index: " +
								// part.matchingMCPartIndex);	
								
								if (part.pid() == LundPID.Pion.lundCode() || part.pid()==LundPID.Kaon.lundCode()) {
									
									for (int j = 0; j < generic_Event.count(); j++) {
										MyParticle part2 = (MyParticle) generic_Event.getParticle(j);
										// Systefm.out.println("lookign at pid " + part2.pid());
										if (part2.pid() == ((-1)*LundPID.Pion.lundCode()) || part2.pid()==((-1)*LundPID.Kaon.lundCode())){
											
											
										
										}
									}
							
								}
							}
							
						}
					}
				}
			}
		}
		this.novel_fitter.saveHitHistograms();
	}
	protected void associateMCWithData(PhysicsEvent generic_Event, PhysicsEvent generic_EventMC) {
		//System.out.println("associate");
				for (int i = 0; i < generic_Event.count(); i++) {
					double minMomDiff = 1000.0;
					int minMomDiffIndex = -1;
					MyParticle part = (MyParticle) generic_Event.getParticle(i);
					double px = part.px();
					double py = part.py();
					double pz = part.pz();
					System.out.println("data part px " + px + " py " + py + " pz " + pz);

					double mom = Math.sqrt(px * px + py * py + pz * pz);
					double theta = Math.toDegrees(Math.atan2(Math.sqrt(px * px + py * py), pz));
					double phi = Math.toDegrees(Math.atan2(py, px));
					System.out.println("data theta "+ theta + " phi " + phi + " mom " +mom);
					for (int j = 0; j < generic_EventMC.count(); j++) {
						MyParticle partMC = (MyParticle) generic_EventMC.getParticle(j);
						double pxMC = partMC.px();
						double pyMC = partMC.py();
						double pzMC = partMC.pz();

						double momMC = Math.sqrt(pxMC * pxMC + pyMC * pyMC + pzMC * pzMC);
						double thetaMC = Math.toDegrees(Math.atan2(Math.sqrt(pxMC * pxMC + pyMC * pyMC), pzMC));
						double phiMC = Math.toDegrees(Math.atan2(pyMC, pxMC));
						//if (i == 0) 
						{
							 System.out.println("MC part px " + pxMC + " py " + pyMC + " pz " + pzMC);
							 System.out.println("mc teta "+ thetaMC + " phi " + phiMC + " mom " +momMC);
						}

						// System.out.println("Looking at MC part "+ j + " rel momDiff: " + ((mom-momMC)/momMC));
					//	if (Math.abs(mom - momMC) < 0.025 * momMC && Math.abs(theta - thetaMC) < 1.0 && Math.abs(phi - phiMC) < 5)
						//let's increase this
						if (Math.abs(mom - momMC) < (0.1 * momMC) && Math.abs(theta - thetaMC) < 2.0 && Math.abs(phi - phiMC) < 10)
						{
							if (Math.abs(mom - momMC) < minMomDiff) {
								minMomDiff = Math.abs(mom - momMC);
								minMomDiffIndex = j;
							}
						}
					}
				System.out.println("associate mc part "+minMomDiffIndex + " with "+ i);
					part.matchingMCPartIndex = minMomDiffIndex;
				}

			}
	
	void doDiHadrons(PhysicsEvent generic_Event, PhysicsEvent generic_EventMC, NovelBaseFitter m_novel_fitter, NovelBaseFitter m_novel_fitterMC) {
		//System.out.println("in do dihad with "+generic_Event.count() + "particles ");
		
		if(novel_fitter.evtNumber==51831184)
		{
			System.out.println("found our run..");
		}
		for (int i = 0; i < generic_Event.count(); i++) 
		
		
		{
			MyParticle part = (MyParticle) generic_Event.getParticle(i);
			int sec= part.FTOFsector;
			if(sec<=0)
				sec=6;
			
			//System.out.println("sec: " + sec+ ", beta is " + part.beta);
			if(part.beta<10.0 && part.beta>-10.0)
			{
				
			}
			//System.out.println("time is: " + part.FTOFTime + " sector: " + part.FTOFsector);
			// System.out.println("matching mc particle index: " +
			// part.matchingMCPartIndex);

			
			 
			
			if (part.pid() == LundPID.Pion.lundCode() || part.pid()==LundPID.Kaon.lundCode()) {
				
				for (int j = 0; j < generic_Event.count(); j++) {
					MyParticle part2 = (MyParticle) generic_Event.getParticle(j);
					// Systefm.out.println("lookign at pid " + part2.pid());
					if (part2.pid() == ((-1)*LundPID.Pion.lundCode()) || part2.pid()==((-1)*LundPID.Kaon.lundCode())){
						
						System.out.println("part mom1: "+ part.p() +" mom2: " +part2.p());
					if(part.e()<1.0 || part2.e()<1.0)
						continue;
					if(part.theta()<0.17 || part2.theta()<0.17)
						continue;
						
					
						
						this.numAllPairs++;
						if (part.matchingMCPartIndex != -1 && part2.matchingMCPartIndex != -1) {
							//System.out.println("found di hadron  candidate with matching MC!!");

							this.numMatch++;
						}
					
					//System.out.println("adding hadron pair");	
					//System.out.println("pair data before: " + this.currentEvent.pairData.size());	
						
					//	System.out.println("pair data now: " + this.currentEvent.pairData.size());	
						
					}
				}
			}
		}

	}
	
}
