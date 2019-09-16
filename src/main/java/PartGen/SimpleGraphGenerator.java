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
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.silent.AtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;


public class SimpleGraphGenerator {
	public static IChemObjectBuilder builder =SilentChemObjectBuilder.getInstance();
	public static IMolecularFormula formula=null;
	public static IAtomContainer acontainer;
	public static Map<Integer, Integer> capacities;
	public static boolean verbose = false;
	public static int atomCount;
	public static int totalCapacity;
	public static int totalValences;
	public static int hydrogens;
	public static int isotopes;
	public static int[] capacity;
	public static int[] valences;
	
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
	 * Add an element ot an int array
	 * @param a int array
	 * @param e int value to add
	 * @return new int array
	 */
	
	public static int[] addElement(int[] a, int e) {
        a  = Arrays.copyOf(a, a.length + 1);
        a[a.length - 1] = e;
        return a;
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
	
	/**
	 * Sum all entires of an int array
	 * @param array int array
	 * @return int the sum
	 */
	public static int sum(int[] array) {
		int sum=0;
		for(int i=0;i<array.length;i++) {
			sum=sum+array[i];
		}
		return sum;
	}
	
	/**
	 * Sum all entires of an int matrix
	 * @param array int array
	 * @return int the sum
	 */
	
	public static int sum(int[][] array) {
		int sum=0;
		for(int i=0;i<array.length;i++) {
			sum=sum+sum(array[i]);
		}
		return sum;
	}
	
	/**
	 * Copy an int matrix
	 * @param matrix int matrix
	 * @return int matrix
	 */
	public static int[][] copy(int[][] matrix){
		int[][] copy = Arrays.stream(matrix).map(r -> r.clone()).toArray(int[][]::new);
		return copy;
	}
	
	/**
	 * Get a column from a matrix
	 * @param matrix an int matrix 
	 * @param index index of a column
	 * @return int array, column
	 */
	
	public static int[] getColumn(int[][] matrix, int index){
	    int size=matrix[0].length;
		int[] column = new int[size]; 
	    for(int i=0; i<size; i++){
	       column[i] = matrix[i][index];
	    }
	    return column;
	}
	
	/**
	 * Adding zeros to the front of an array
	 * @param array int array
	 * @param zeros number of zeros to add
	 * @return updated int array
	 */
	
	public static int[] addZerosF(int[] array, int zeros) {
		int[] arr= new int[zeros];
		for(int i=0;i<zeros;i++) {
			arr[i]=0;
		}
		for(int i=0;i<array.length;i++) {
			arr=addElement(arr, array[i]);
		}
		return arr;
	}
	
	/**
	 * Summing valences of all the atoms
	 * @param ac atomcontainer
	 * @return int summation of all valences.
	 */
	
	public static int totalValences(IAtomContainer ac) {
		int val=0;
		for(IAtom atom: acontainer.atoms()) {
			val+=capacities.get(atom.getAtomicNumber())+1;
		}
		return val;
	}
	
	/**
	 * Valence check for an adjacecny matrix
	 * @param mat int matrix
	 * @param valence valence of an atom
	 * @return boolean 
	 */
	
	public static boolean valenceCheck(int[][] mat){
		boolean check=true;
		for(int i=0;i<mat.length;i++) {
			if(!(sum(mat[i])<=valences[i] && sum(getColumn(mat,i))<=valences[i])) {
				check=false;
				break;
			}
		}
		return check;	
	}
	
	/**
	 * Setting general parameters of the class.
	 * @param formula Molecular formula.
	 * @return
	 */
	
	public static int[] setValues(IMolecularFormula formula) {
		int size= formula.getIsotopeCount();
		int[] capacity = new int[size];
		int[] valences = new int[size];
		int[] totalAtom = new int[size];
		int i=0;
		for(IIsotope top:formula.isotopes()) {
			Integer atomNo=top.getAtomicNumber();
			totalAtom[i]=formula.getIsotopeCount(top);
			valences[i]=capacities.get(atomNo)+1;
			capacity[i]=capacities.get(atomNo)*formula.getIsotopeCount(top);
			i++;
		}
		
		SimpleGraphGenerator.capacity=capacity;
		SimpleGraphGenerator.valences=valences;
		SimpleGraphGenerator.totalValences=totalValences(acontainer);
		return capacity;
	}
	
	/**
	 * Integer partitioning
	 * @param index row number in the matrix
	 * @param n integer to distribute
	 * @param fill sum of filled entries
	 * @param d 
	 * @param depth
	 * @return
	 */
	
	public static List<int[]> partition(int index, int n,int fill, int d,int depth) {
		if(d==depth) {
			List<int[]> array= new ArrayList<int[]>();
			int[] take=new int[0];
			array.add(take);
			return array;
		}
		return buildArray(index, n,fill,d,depth);	
	}
	
	public static List<int[]> buildArray(int index, int n,int fill,int d, int depth){
		List<int[]> array= new ArrayList<int[]>();
		for(int i=0;i<2;i++) {
			for(int[] item: partition(index, n-i,fill,d,depth+1)) {
				if(item.length==0) {
					item=addElement(item,i);
					if(sum(item)<=valences[item.length-1]) {
						array.add(item);
					}
				}else {
					if(item[item.length-1]<=i) {
						item=addElement(item,i);
						if(sum(item)<=valences[item.length-1]) {
							array.add(item);
						}
					}
				}
			}
		}
		return array;
	}
	
	/**
	 * Generating adjacency matrices of an atomcontainer
	 * @param ac atomcontainer
	 * @param index row index of adj. matrix
	 * @param bonds number of bonds to distribute. 
	 * @param matrices list of filled matrices.
	 * @param output the list of fully filled adj. matrices
	 * @return
	 */
	
	public static List<int[][]> generate(IAtomContainer ac,int index,int bonds,List<int[][]> matrices,List<int[][]> output) {
		int atomSize= ac.getAtomCount();
		if(index==atomSize-1) {
			for(int[][] mat:matrices) {
				if(valenceCheck(mat)) {
					System.out.println(Arrays.deepToString(mat));
					output.add(mat);
				}
			}
		}else {
			for(int[][] mat:matrices) {
				int de=atomSize-(index+1);
				if(de!=0) {
					for(int[] array:partition(index,atomCount-index,sum(mat[index]),de,0)){
						if(sum(array)!=0 && sum(array)<=valences[index]) {
							List<int[][]> list= new ArrayList<int[][]>();
							int zeros=atomSize-array.length;
							array=addZerosF(array,zeros);
							int[][] copy=copy(mat);
							copy[index]=array;
							list.add(copy);
							generate(ac,index+1,bonds-sum(copy),list,output);
						}
					}
				}
			}
		}
		return matrices;
	}
	
	/**
	 * Generate atomcontainers for a list of adjacency matrices.
	 * @param adcacencyMatrices 
	 * @return list of atomcontainers
	 * @throws CloneNotSupportedException
	 */
	
	public static List<IAtomContainer> generateAtomContainers(List<int[][]> adcacencyMatrices) throws CloneNotSupportedException{
		List<IAtomContainer> acontainers= new ArrayList<IAtomContainer>();
		for(int[][] adjacency:adcacencyMatrices) {
			IAtomContainer ac=acontainer.clone();
			acontainers.add(setBonds(ac,adjacency));
		}
		return acontainers;
	}
	
	/**
	 * Setting bonds to an atomContainer
	 * @param ac
	 * @param adjacency
	 * @return
	 * @throws CloneNotSupportedException
	 */
	
	public static IAtomContainer setBonds(IAtomContainer ac,int[][] adjacency) throws CloneNotSupportedException {
		IAtomContainer ac2=ac.clone();
		for(int i=0;i<adjacency.length;i++) {
			for(int j=i+1;j<adjacency.length;j++) {
				ac2.addBond(i, j, Order.SINGLE);
			}
		}
		return ac2;
	}
	
	/**
	 * To initialise the inputs and run the functions while recording the duration time.
	 * @throws CDKException 
	 */
	
	public static void run(IMolecularFormula formula) throws FileNotFoundException, UnsupportedEncodingException, CloneNotSupportedException, CDKException {
		long startTime = System.nanoTime(); //Recording the duration time.
		int hydrogen=formula.getIsotopeCount(builder.newInstance(IIsotope.class, "H"));
		String formulaString =MolecularFormulaManipulator.getString(formula);
		SimpleGraphGenerator.isotopes=formula.getIsotopeCount()-1;
		formula.removeIsotope(builder.newInstance(IIsotope.class, "H"));
		IAtomContainer ac=MolecularFormulaManipulator.getAtomContainer(formula);
		SimpleGraphGenerator.atomCount=ac.getAtomCount();
		SimpleGraphGenerator.acontainer=ac;
		SimpleGraphGenerator.hydrogens=hydrogen;
		setValues(formula);	
		if(verbose) {
			System.out.println("For molecular formula "+ formulaString +", generating structures...");
		}
		List<int[][]> matrices= new ArrayList<int[][]>();
		for(int[] array:partition(0,atomCount-1,0,atomCount-1,0)){
			if(sum(array)!=0 && sum(array)<=valences[0]) { 
				//TODO: first line can be zero 
				int zeros=atomCount-array.length;
				array=addZerosF(array,zeros);
				int[][] mat= new int[atomCount][atomCount];
				mat[0]=array;
				matrices.add(mat);
			}
		}

		List<int[][]> output= new ArrayList<int[][]>();
		generate(ac,1,(totalValences-hydrogens)/2,matrices,output);
		for(int[][] mat: output) {
			System.out.println(Arrays.deepToString(mat));
		}
		if(verbose) {
			System.out.println("Number of structures: "+output.size());
		}
		long endTime = System.nanoTime()- startTime;
        double seconds = (double) endTime / 1000000000.0;
		DecimalFormat d = new DecimalFormat(".###");
		if(verbose) {
			System.out.println("Duration:"+" "+d.format(seconds));
		}
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
			String header = "\nFor a molecular formula, it generates all possible structures";
			String footer = "\nPlease report issues at https://github.com/MehmetAzizYirik/PartGen";
			formatter.printHelp( "java -jar partgen.jar", header, options, footer, true );
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
	
	public static void main(String[] arguments) throws FileNotFoundException, UnsupportedEncodingException, CloneNotSupportedException {
		Generator gen= new Generator();
		String[] arguments1= {"-f","C6H6","-v"};
		try {
			gen.parseArguments(arguments1);
			Generator.run(Generator.formula);
		} catch (Exception e) {
			if (Generator.verbose) e.getCause(); 
		}
	}
}
