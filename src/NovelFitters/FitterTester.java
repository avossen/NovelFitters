package NovelFitters;

import java.io.File;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.nio.file.PathMatcher;
import java.nio.file.Paths;

import org.jlab.clas.physics.EventFilter;
import org.jlab.clas.physics.PhysicsEvent;
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
	}

	
	public void analyze(String[] args) {
		
		NovelBaseFitter.useStefanElectronCuts=true;
		NovelBaseFitter.useStefanHadronCuts=true;
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
						//System.out.println("new event----\n\n");
						HipoDataEvent event = (HipoDataEvent) reader.getNextEvent();
						
						// apply fitter to your data event
						PhysicsEvent generic_Event = novel_fitter.getPhysicsEvent(event);
						PhysicsEvent generic_EventMC=new PhysicsEvent();
						if(isMC)
						{
							generic_EventMC = novel_fitterMC.getPhysicsEvent(event);
						}
						
						// novel_fitter.Walt);
						if(generic_Event.count()>2)
							System.out.println("looking at event with " + generic_Event.count() + " particles ");
						if (filter.isValid(generic_Event) == true) { // apply filter to current event
							// look at all particles
							System.out.println("we have " +generic_Event.count());
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
	}
	
}
