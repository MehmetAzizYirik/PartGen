package PartGen;
import java.io.FileNotFoundException;
import java.io.UnsupportedEncodingException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.IntStream;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IIsotope;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.silent.AtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;


public class Generator {
	public static IChemObjectBuilder builder =SilentChemObjectBuilder.getInstance();
	public static IMolecularFormula formula=null;
	public static IAtomContainer acontainer;
	public static Map<Integer, Integer> capacities;
	public static boolean verbose = false;
	public static int size;
	public static int totalCapacity;
	public static int totalValences;
	public static int isotopes;
	public static int[] capacity;
	public static int[] valences;
	public static int[][] maxMultiplicity;
	public static int totalHydrogen; // Total number of hydrogens.
	public static int[] totalAtom; // Total number of atoms.
	public static int hydrogens2distribute;
	
	static {
		//The atom capacities from MOLGEN book. Capacity of an atom equals to 
		capacities = new HashMap<Integer, Integer>();
		capacities.put(6, 3);
		capacities.put(7, 2);
		capacities.put(8, 1);
		capacities.put(16, 1);
		capacities.put(15, 2);
		capacities.put(9, 0);
		capacities.put(53, 0);
		capacities.put(1, 0);
		capacities.put(17, 0);
		capacities.put(35, 0);
		capacities.put(53, 0);
		
	}
	
	/**
	 * The basic functions used in the hydrogen distributor.
	 */
	
	public static int[] addElement(int[] a, int e) {
        a  = Arrays.copyOf(a, a.length + 1);
        a[a.length - 1] = e;
        return a;
    }
	
	public static int[] setValues(IMolecularFormula formula) {
		int size= formula.getIsotopeCount();
		int[] capacity = new int[size];
		int[] valences = new int[size];
		int[] totalAtom = new int[size];
		int i=0;
		for(IIsotope top:formula.isotopes()) {
			Integer atomNo=top.getAtomicNumber();
			totalAtom[i]=formula.getIsotopeCount(top);
			valences[i]=capacities.get(atomNo);
			capacity[i]=capacities.get(atomNo)*formula.getIsotopeCount(top);
			i++;
		}
		
		Generator.capacity=capacity;
		Generator.valences=valences;
		Generator.totalCapacity=sum(capacity);
		Generator.totalValences=totalValences(acontainer);
		Generator.totalAtom=totalAtom;
		Generator.maxMultiplicity= maxMultiplicity(valences);
		return capacity;
	}
	
	public static int totalValences(IAtomContainer ac) {
		int val=0;
		for(IAtom atom: acontainer.atoms()) {
			val+=capacities.get(atom.getAtomicNumber())+1;
		}
		return val;
	}
	
	public static int[][] maxMultiplicity(IAtomContainer ac) {
		int atoms= ac.getAtomCount();
		int comb= combination(atoms,2);
		int[][] mult= new int[atoms][atoms];
		int count=0;
		for(int i=0;i<atoms;i++) {
			for(int j=i+1;j<atoms;j++) {
				if(ac.getAtom(i).getSymbol()!=ac.getAtom(j).getSymbol()) {
					mult[i][j]=Math.min(capacities.get(ac.getAtom(i).getAtomicNumber()), capacities.get(ac.getAtom(j).getAtomicNumber()));
				}else {
					mult[i][j]=Math.min(capacities.get(ac.getAtom(i).getAtomicNumber()), capacities.get(ac.getAtom(j).getAtomicNumber()))-1;
				}
				count++;
			}
		}
		return mult;
	}
	
	public static int factorial(int i){
		if (i==0){
			return 1;
		}
		return i * factorial(i - 1);
	}
	
	public static int combination(int m, int n) {
		return factorial(m) / (factorial(n) * factorial(m - n));
	}
	
	public static int sum(int[] array) {
		int sum=0;
		for(int i=0;i<array.length;i++) {
			sum=sum+array[i];
		}
		return sum;
	}
	
	public static int sum(int[][] array) {
		int sum=0;
		for(int i=0;i<array.length;i++) {
			sum=sum+sum(array[i]);
		}
		return sum;
	}
	
	public static int[][] copy(int[][] matrix){
		int[][] copy = Arrays.stream(matrix).map(r -> r.clone()).toArray(int[][]::new);
		return copy;
	}
	
	public static int[] mergeArrays(List<int[]> arrays) {
		int size = 0;
		for (int[] array : arrays) {
			size += array.length;
		}
		int[] mergedArray = new int[size];
		int index = 0;
		for (int[] array : arrays) {
			for (int i : array) {
				mergedArray[index++] = i;
		    }
		}
		return mergedArray;
	}
	
	public static int[] arraySum(int[] a, int[] b) {
		List<int[]> arrays= new ArrayList<int[]>();
		arrays.add(a);
		arrays.add(b);
		return mergeArrays(arrays);
	}
	
	public static List<List<int[]>> buildLists(int n){
		List<List<int[]>> lists= new ArrayList<List<int[]>>();
		for (int i=0; i<n; ++i) {
			List<int[]> ilist= new ArrayList<int[]>();
			lists.add(ilist);
		}
		return lists;
	}
	public static Set<int[][]> matr= new HashSet<int[][]>();
	public static List<int[][]> generate(IAtomContainer ac,int index,int bonds,List<int[][]> matrices) {
		int atomSize= ac.getAtomCount();
		if(index==atomSize-1) {
			for(int[][] mat:matrices) {
				matr.add(mat);
				//System.out.println(Arrays.deepToString(mat));
			}
			System.out.println(matr.size());
		}else {
			List<int[][]> list= new ArrayList<int[][]>();
			int capa=capacities.get(ac.getAtom(index).getAtomicNumber());
			for(int[][] mat:matrices) {
				int dis=Math.min(capa, bonds-sum(mat));
				int de=atomSize-(index+1);
				if(de!=0) {
					for(int[] dene:partition(index,dis,de,0)){
						if(sum(dene)!=0) {
							int[][] copy=copy(mat);
							copy[index]=dene;
							list.add(copy);
						}
					}
				}
			}
			generate(ac,index+1,bonds,list);
		}
		return matrices;
	}
	/**
	 * To initialise the inputs and run the functions while recording the duration time.
	 * @throws CDKException 
	 */
	
	public static void run(IMolecularFormula formula) throws FileNotFoundException, UnsupportedEncodingException, CloneNotSupportedException, CDKException {
		long startTime = System.nanoTime(); //Recording the duration time.
		int hydrogen=formula.getIsotopeCount(builder.newInstance(IIsotope.class, "H"));
		String formulaString =MolecularFormulaManipulator.getString(formula);
		Generator.isotopes=formula.getIsotopeCount()-1;
		formula.removeIsotope(builder.newInstance(IIsotope.class, "H"));
		IAtomContainer ac=MolecularFormulaManipulator.getAtomContainer(formula);
		Generator.size=ac.getAtomCount();
		Generator.acontainer=ac;
		setValues(formula);
		Generator.totalHydrogen=hydrogen;		
		if(verbose) {
			System.out.println("For molecular formula "+ formulaString +", calculating all the possible distributions of "+totalHydrogen+" "+"hydrogens ..." );
		}
		List<int[]> result= new ArrayList<int[]>();
		int count=0;
		List<int[]> iarrays= new ArrayList<int[]>();
		int[] array = new int[0];
		distribute(iarrays,(totalValences-totalHydrogen),array);
		count=iarrays.size();
		result= iarrays;
		if(verbose) {
			System.out.println("Number of distributions: "+count);
		}
		long endTime = System.nanoTime()- startTime;
        double seconds = (double) endTime / 1000000000.0;
		DecimalFormat d = new DecimalFormat(".###");
		if(verbose) {
			System.out.println("Duration:"+" "+d.format(seconds));
		}
	}
	
	/**
	 * These functions are built for the integer partitioning problem.
	 */
	
	public static List<int[]> partition(int index, int n, int d,int depth) {
		if(d==depth) {
			List<int[]> array= new ArrayList<int[]>();
			int[] take=new int[0];
			array.add(take);
			return array;
		}
		return buildArray(index, n,d,depth);	
	}
	
	public static int index(int ind) {
		if(ind==0) {
			return 0;
		}else {
			return ind-1;
		}
	}
	
	public static List<int[]> buildArray(int index, int n,int d, int depth){
		List<int[]> array= new ArrayList<int[]>();
		IntStream range = IntStream.rangeClosed(0,n);
		for(int i:range.toArray()) {
			for(int[] item: partition(index, n-i,d,depth+1)) {
				//System.out.println(Arrays.toString(item));
				if(i<=2) { //maxMultipliciy almadi o matrix den.
					item=addElement(item,i);
			        //if(item.length==d) {
			        	//if(sum(item)==valence) {
			        		array.add(item);
			        	//}
			        //}else {
			        	//array.add(item);
			        //}
				}
			}
		}
		return array;
	}
	
	public static int[] addZeros(int[] array, int zeros) {
		for(int i=0;i<zeros;i++) {
			array=addElement(array,0);
		}
		return array;
	}
	
	public static void distribute(List<int[]> arrays,int valence,int[]arr,int atoms) throws CloneNotSupportedException {
		if(valence==0){
			System.out.println("n"+" "+Arrays.toString(arr));
			arrays.add(arr);
		}else { 
			System.out.println(valence+" "+sum(arr));
			int dis= Math.abs(valence-sum(arr));
			System.out.println(dis+" "+Arrays.toString(arr)+" "+"ickisim");
			for(int i = dis; i > 0; i--) {
				distribute(arrays,(valence-i),addElement(arr,i),atoms); 
			}
		}
	}
	
	
	/**
	 * Functions seting the implicit hydrogens of the atoms.
	 * @throws CloneNotSupportedException 
	 */
	
	public static List<IAtomContainer> generateAtomContainers(List<int[]> distributions) throws CloneNotSupportedException{
		List<IAtomContainer> acontainers= new ArrayList<IAtomContainer>();
		for(int[] array:distributions) {
			IAtomContainer ac=acontainer.clone();
			acontainers.add(setHydrogens(ac,array));
		}
		return acontainers;
	}
	
	public static IAtomContainer setHydrogens(IAtomContainer ac,int[] distribution) throws CloneNotSupportedException {
		IAtomContainer ac2=ac.clone();
		for(int i=0;i<distribution.length;i++) {
			ac2.getAtom(i).setImplicitHydrogenCount(distribution[i]);
		}
		return ac2;
	}
	
	void parseArguments(String[] arguments) throws ParseException
	{
		Options options = setupOptions(arguments);	
		CommandLineParser parser = new DefaultParser();
		try {
			CommandLine cmd = parser.parse(options, arguments);
			String formula = cmd.getOptionValue("formula");
			Generator.formula=MolecularFormulaManipulator.getMolecularFormula(formula, builder);
			if (cmd.hasOption("verbose")) Generator.verbose = true;
		
		} catch (ParseException e) {
			HelpFormatter formatter = new HelpFormatter();
			formatter.setOptionComparator(null);
			String header = "\nFor a molecular formula, it calculates all the possible hydrogen distributions to the atoms.";
			String footer = "\nPlease report issues at https://github.com/MehmetAzizYirik/Generator";
			formatter.printHelp( "java -jar Generator.jar", header, options, footer, true );
			throw new ParseException("Problem parsing command line");
		}
	}
	
	private Options setupOptions(String[] arguments)
	{
		Options options = new Options();
		Option formula = Option.builder("f")
			     .required(true)
			     .hasArg()
			     .longOpt("formula")
			     .desc("Molecular Formula (required)")
			     .build();
		options.addOption(formula);	
		Option verbose = Option.builder("v")
			     .required(false)
			     .longOpt("verbose")
			     .desc("Print messages about the distributor")
			     .build();
		options.addOption(verbose);	
		return options;
	}
	public static int total=6;
	public static void main(String[] arguments) throws FileNotFoundException, UnsupportedEncodingException, CloneNotSupportedException {
		/**Generator distribution= new Generator();
		String[] arguments1= {"-f","C6H12","-v"};
		try {
			distribution.parseArguments(arguments1);
			Generator.run(Generator.formula);
		} catch (Exception e) {
			if (Generator.verbose) e.getCause(); 
		}**/
		IMolecularFormula formula= MolecularFormulaManipulator.getMolecularFormula("C6H12", builder);
		int hydrogen=formula.getIsotopeCount(builder.newInstance(IIsotope.class, "H"));
		String formulaString =MolecularFormulaManipulator.getString(formula);
		Generator.isotopes=formula.getIsotopeCount()-1;
		formula.removeIsotope(builder.newInstance(IIsotope.class, "H"));
		IAtomContainer ac=MolecularFormulaManipulator.getAtomContainer(formula);
		Generator.maxMultiplicity=maxMultiplicity(ac);
		int say= ac.getAtomCount();
		int comb= combination(say,2);
		/**List<int[][]> matrices= new ArrayList<int[][]>();
		System.out.println(capacities.get(ac.getAtom(0).getAtomicNumber()));
		for(int[] dene:partition(0,capacities.get(ac.getAtom(0).getAtomicNumber())+1,capacities.get(ac.getAtom(0).getAtomicNumber()),say,0)){
			int[][] mat= new int[say][say];
			mat[0]=dene;
			matrices.add(mat);
		}
		System.out.println(matrices.size());
		for(int a=1;a<ac.getAtomCount();a++) {
			int cap=capacities.get(ac.getAtom(a).getAtomicNumber());
			List<int[][]> matrices2= new ArrayList<int[][]>();
			for(int[][] mat:matrices) {
				for(int[] dene:partition(a,cap+1,cap,say,0)){
					mat[a]=dene;
					matrices2.add(mat);
				}
			}
			System.out.println(matrices2.size());
			matrices=matrices2;
		}**/
		List<int[][]> matrices= new ArrayList<int[][]>();
		for(int[] dene:partition(0,capacities.get(ac.getAtom(0).getAtomicNumber()),say-1,0)){
			if(sum(dene)!=0) {
				//System.out.println(Arrays.toString(dene));
				int[][] mat= new int[say][say];
				mat[0]=dene;
				matrices.add(mat);
			}
		}
		
		generate(ac,1,comb,matrices);
		
	}
}
