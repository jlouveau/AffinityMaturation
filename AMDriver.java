
package tip.maam3;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.apache.commons.collections.primitives.ArrayDoubleList;
import org.apache.commons.collections.primitives.DoubleList;

import tip.app.TipLog;
import tip.app.TipProperties;
import tip.app.TipRuntimeException;
import tip.io.FileUtil;
import tip.math.MeanEstimator;
import tip.math.StatUtil;
import tip.math.TipRandom;
import tip.math.TipVector;
import tip.math.BitVector;
import tip.math.DoubleUtil;

public final class AMDriver {
	private final File propFile;
	private final File analysisFolder;
	private final File rootDir;

	private final int randomSeed;
	private final int trialLimit;
	private final int plasmaTarget;

	private final AntigenProp agProp;
	private final BCellProp bcProp;
	private final GerminalCenterProp gcProp;
	private final VaccinationStrategyProp vsProp;
	private final ViralChallengeProp vcProp;

	// Unique plasma cells produced by affinity maturation...
	private final Collection<BCell> plasmaCells = new HashSet<BCell>();

	// Nearest epitope for each plasma cell...
	private final Map<BCell, Epitope> targetEpitopes = new HashMap<BCell, Epitope>();

	// Viral challenge to each plasma cell...
	private final Map<BCell, ViralChallenge> viralChallenges = new HashMap<BCell, ViralChallenge>();

	// Production rate for each trial...
	private final ArrayDoubleList productionRates = new ArrayDoubleList();

	// Current trial information...
	private int trialIndex = 0;
	private GerminalCenter activeGC;
	private int extinguishedGC = 0;
	private double survivalRate;

	private PrintWriter trialDetailWriter;

	// Germinal centers produced for each trial in the simulation...
	// private final List<GerminalCenter> germinalCenters = new
	// ArrayList<GerminalCenter>();

//	private static String nameDirectory;

	private static final String TRIAL_SUMMARY_BASENAME = "trial-summary.csv";
	private static final String PLASMA_DETAIL_BASENAME = "plasma-detail.csv";
	private static final String PLASMA_SUMMARY_BASENAME = "plasma-summary.csv";

	private AMDriver(File propFile, File analysisFolder) {
		this.propFile = propFile;
		this.analysisFolder = analysisFolder;
		this.rootDir = FileUtil.getParentFile(propFile);
		
		//this.nameDirectory = null;
		
		TipProperties.loadFile(propFile);

		this.randomSeed = TipProperties.getRequiredInt("AMDriver.randomSeed");
		this.trialLimit = TipProperties.getRequiredInt("AMDriver.trialLimit");
		this.plasmaTarget = TipProperties.getRequiredInt("AMDriver.plasmaTarget");

		this.agProp = AntigenProp.instance();
		this.bcProp = BCellProp.instance();
		this.gcProp = GerminalCenterProp.instance();
		this.vcProp = ViralChallengeProp.instance();
		this.vsProp = VaccinationStrategyProp.instance();
	}

	/**Creates AMDriver with propFile and analysisFolder.
	 * 
	 * @param propFile is the property file.
	 * 
	 * @param analysisFolder is the folder in which the many trials for a given propFile are analyzed.
	 * it is named after the property file + index if a simulation for that property file has already been analyzed
	 */
	
	public static void run(File propFile, File analysisFolder) {
		AMDriver driver = new AMDriver(propFile, analysisFolder);
		driver.run();
	}

	private void run() {
		TipRandom.initialize(randomSeed);

		while (continueTrials()){
			runTrial();
		}
		trialDetailWriter.close();
		
		survivalRate = 1 - DoubleUtil.ratio(extinguishedGC, trialIndex);
		
		if (!plasmaCells.isEmpty())
			analyzeTrials();

	}
	
	private boolean continueTrials() {
		return trialIndex < trialLimit && plasmaCells.size() < plasmaTarget;
	}

	private void runTrial() {
		TipLog.info("TRIAL: %4d", trialIndex);
				
		trialDetailWriter =	createWriter("trialIndex" + String.valueOf(trialIndex) +".csv");
		writeTrialDetailHeader();
		
		activeGC = GerminalCenter.run(trialDetailWriter);
		trialDetailWriter.flush();
		
		if (activeGC.IsExtinguished()) {
			extinguishedGC++;
		}

		plasmaCells.addAll(activeGC.getPlasmaCells());
		targetEpitopes.putAll(activeGC.mapTargets());
		viralChallenges.putAll(activeGC.mapChallenges());

		productionRates.add(activeGC.computePlasmaProductionRate());

		++trialIndex;
	}

	private double getSurvivalRate() {
		return this.survivalRate;
	}

	private void analyzeTrials() {

		writePlasmaDetail();
		writePlasmaSummary();
	}

	private void writeTrialDetailHeader() {
		trialDetailWriter.print("GCcycleIndex");
		trialDetailWriter.print(",");
		trialDetailWriter.print("#recycled B cells");
		trialDetailWriter.print(",");
		trialDetailWriter.print("Mean #deleterious mutations");
		trialDetailWriter.print(",");
		trialDetailWriter.print("Mean #beneficial mutations");
		trialDetailWriter.print(",");
		trialDetailWriter.print("Mean #total mutations");

		trialDetailWriter.print(",");
		trialDetailWriter.print("Vaccine concentration");
		trialDetailWriter.print(",");
		trialDetailWriter.print("Mean max energy change from founder");
		trialDetailWriter.print(",");
		trialDetailWriter.print("Variance max energy change");
		trialDetailWriter.print(",");
		trialDetailWriter.print("Total Affinity change");

		trialDetailWriter.println();
	}

	private void writePlasmaDetail() {
		PrintWriter writer = createWriter(PLASMA_DETAIL_BASENAME);
		writePlasmaDetailHeader(writer);

		for (BCell plasmaCell : plasmaCells)
			writePlasmaDetail(writer, plasmaCell);

		writer.close();
	}

	private void writePlasmaDetailHeader(PrintWriter writer) {
		// int ConservedLength = agProp.getConservedLength();
		// int VariableLength = agProp.getVariableLength();
		// int epitopeLength = getEpitopeLength();

		writer.print("trialIndex");
		writer.print(",");
		writer.print("generation");
		writer.print(",");
		writer.print("mutationCount");
		writer.print(",");
		writer.print("ConservedLength");
		writer.print(",");
		writer.print("VariableLength");

		writer.print(",");
		writer.print("threshold");
		writer.print(",");
		writer.print("breadth");

		writer.print(",");
		writer.print("targetAffinity");
//		writer.print(",");
//		writer.print("targetSpecificity");

		writer.println();
	}

	private void writePlasmaDetail(PrintWriter writer, BCell plasmaCell) {
		int conservedLength = agProp.getConservedLength();
		int variableLength = agProp.getVariableLength();


		Epitope targetEpitope = targetEpitopes.get(plasmaCell);
		DecimalFormat formatter = new DecimalFormat("##0.0#####");
		ViralChallenge challenge = viralChallenges.get(plasmaCell);

		trialDetailWriter.print(trialIndex);
		trialDetailWriter.print(",");
		writer.print(plasmaCell.getGeneration());
		writer.print(",");
		writer.print(plasmaCell.getMutationCount());
		writer.print(",");
		writer.print(conservedLength);
		writer.print(",");
		writer.print(variableLength);

			writer.print(",");
			writer.print(formatter.format(ViralChallenge.getThreshold()));
			writer.print(",");
			writer.print(formatter.format(challenge.getBreadth()));


		writer.print(",");
		writer.print(formatter.format(plasmaCell.getReceptor().calculateEnergy(targetEpitope)));
		// writer.print(",");
		// writer.print(formatter.format(Receptor.computeSpecificity(plasmaCell.getReceptor(),targetEpitope)));

		writer.println();
	}

	// private int getEpitopeLength() {
	// return agProp.getConservedLength() + agProp.getVariableLength();
	// }

	private void writePlasmaSummary() {
		PrintWriter writer = createWriter(PLASMA_SUMMARY_BASENAME);
		DecimalFormat formatter = new DecimalFormat("##0.0#####");

		writer.print("vaccineCount");
		writer.print(",");
		writer.print("conservedLength");
		writer.print(",");
		writer.print("variableLength");
		writer.print(",");
		writer.print("generationMean");
		writer.print(",");
		writer.print("mutationCountMean");
		writer.print(",");
		writer.print("Survival Rate");

			writer.print(",");
			writer.print("neutThreshold");
			writer.print(",");
			writer.print("breadthMean");

		writer.println();

		writer.print(vsProp.getVaccineCount());
		writer.print(",");
		writer.print(agProp.getConservedLength());
		writer.print(",");
		writer.print(agProp.getVariableLength());
		writer.print(",");
		writer.print(formatter.format(computeGenerationMean(plasmaCells)));
		writer.print(",");
		writer.print(formatter.format(computeMutationCountMean(plasmaCells)));
		writer.print(",");
		writer.print(formatter.format(getSurvivalRate()));

			writer.print(",");
			writer.print(vcProp.getAffinityThreshold());
			writer.print(",");
			writer.print(formatter.format(computeBreadthMean(plasmaCells)));

		writer.println();
		writer.close();
	}

	private static double computeGenerationMean(Collection<BCell> cells) {
		ArrayDoubleList generations = new ArrayDoubleList();

		for (BCell cell : cells)
			generations.add(cell.getGeneration());

		return StatUtil.mean(generations.toArray());
	}

	private static double computeMutationCountMean(Collection<BCell> cells) {
		ArrayDoubleList mutationCounts = new ArrayDoubleList();

		for (BCell cell : cells)
			mutationCounts.add(cell.getMutationCount());

		return StatUtil.mean(mutationCounts.toArray());
	}

	private double computeBreadthMean(Collection<BCell> cells) {
		ArrayDoubleList breadths = new ArrayDoubleList();

		for (BCell cell : cells)
			breadths.add(viralChallenges.get(cell).getBreadth());

		return StatUtil.mean(breadths.toArray());
	}

	private PrintWriter createWriter(String baseName) {
		//
		// Rather not pollute the code with "throws IOException"
		// everywhere we use a PrintWriter...
		//
		File file = new File(analysisFolder, baseName);
		System.out.println("create writer in folder named " + analysisFolder);
		try {
			TipLog.info("Writing [%s]...", file);
			return new PrintWriter(new FileWriter(file));
		} catch (IOException ex) {
			throw new TipRuntimeException(ex);
		}
	}

	private static void usage() {
		System.err.println("Usage: java [...] configDir and analysisDir");
		System.exit(1);
	}

	public static void main(String[] args) {
		if (args.length != 2)
			usage();

		String configDirName = args[0];
		String analysisDirName = args[1];
		File config = null; 
		File[] paths;

		// create new files
		config = new File(configDirName);
		if (config.exists()) {

			// returns pathnames for files and directory
			paths = config.listFiles();

			// for each pathname in pathname array
			//for(File propFile:paths)
			File propFile = new File("/Users/jlouveau/Dropbox (MIT)/AKC_group/Config_AM/SingleAntigen/test.prop");
			{
				
				// prints file and directory paths
				System.out.println("running configuration file " + propFile);	
				String property = propFile.getName();
				String analysisFullPath = "";
				String analysisMainDir = analysisDirName + "/" + property.substring(property.lastIndexOf("/") + 1);
				File analysisFolder = new File(analysisMainDir);
				int i = 0;
				while (analysisFolder.exists()){
					i++;
					analysisFullPath = analysisMainDir + Integer.toString(i);
					analysisFolder = new File(analysisFullPath);
				}

				boolean bool = false;
				// create analysis directory
				bool = analysisFolder.mkdir();
				if ( !bool ){
					throw new IllegalStateException("Cannot create analysis directory");
				}

				run(propFile, analysisFolder);
			}
		}
		else{
			throw new IllegalStateException("No configuration directory");
		}        
	}
}
