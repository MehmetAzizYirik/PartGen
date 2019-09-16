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
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IIsotope;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.silent.AtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;


public class HydrocarbonGenerator {
	public static IChemObjectBuilder builder =SilentChemObjectBuilder.getInstance();
	public static IMolecularFormula formula=null;
	public static IAtomContainer acontainer;
	public static boolean verbose = false;
	public static int atomCount;
	public static int bonds;
	public static int isotopes;
	public static int hydrogens; 
	
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
	 * Valence check for an adjacecny matrix
	 * @param mat int matrix
	 * @param valence valence of an atom
	 * @return boolean 
	 */
	
	public static boolean valenceCheck(int[][] mat, int valence){
		boolean check=true;
		for(int i=0;i<mat.length;i++) {
			if(!(sum(mat[i])<=valence && sum(getColumn(mat,i))<=valence)) {
				check=false;
				break;
			}
		}
		return check;	
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
		for(int i=Math.min(4-fill,3);i>=0;i--) {
			for(int[] item: partition(index, n-i,fill,d,depth+1)) {
				if(item.length==0) {
					item=addElement(item,i);
					if(sum(item)<=4) {
						array.add(item);
					}
				}else {
					if(item[item.length-1]<=i) {
						item=addElement(item,i);
						if(sum(item)<=4) {
							array.add(item);
						}
					}
				}
			}
		}
		return array;
	}
	public static List<int[][]> generate(IAtomContainer ac,int index,int bonds,List<int[][]> matrices,List<int[][]> output) {
		int atomSize= ac.getAtomCount();
		if(index==atomSize-1) {
			for(int[][] mat:matrices) {
				if(valenceCheck(mat,4)) {
					output.add(mat);
				}
			}
		}else {
			for(int[][] mat:matrices) {
				int dis=bonds-sum(mat);
				int de=atomSize-(index+1);
				if(de!=0) {
					for(int[] dene:partition(index,dis,sum(mat[index]),de,0)){
						if(sum(dene)!=0 && sum(dene)<=4) {
							List<int[][]> list= new ArrayList<int[][]>();
							int zeros=atomSize-dene.length;
							dene=addZerosF(dene,zeros);
							int[][] copy=copy(mat);
							copy[index]=dene;
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
	 * To initialise the inputs and run the functions while recording the duration time.
	 * @throws CDKException 
	 */
	
	public static void run(IMolecularFormula formula) throws FileNotFoundException, UnsupportedEncodingException, CloneNotSupportedException, CDKException {
		long startTime = System.nanoTime(); //Recording the duration time.
		int hydrogen=formula.getIsotopeCount(builder.newInstance(IIsotope.class, "H"));
		String formulaString =MolecularFormulaManipulator.getString(formula);
		HydrocarbonGenerator.isotopes=formula.getIsotopeCount()-1;
		HydrocarbonGenerator.hydrogens=hydrogen;
		formula.removeIsotope(builder.newInstance(IIsotope.class, "H"));
		IAtomContainer ac=MolecularFormulaManipulator.getAtomContainer(formula);
		HydrocarbonGenerator.atomCount=ac.getAtomCount();
		HydrocarbonGenerator.bonds= atomCount-1;
		HydrocarbonGenerator.acontainer=ac;		
		if(verbose) {
			System.out.println("For molecular formula "+ formulaString +", generating hydrocarbons...");
		}
		List<int[][]> matrices= new ArrayList<int[][]>();
		for(int[] array:partition(0,4,0,atomCount-1,0)){
			if(sum(array)!=0 && sum(array)<=4) {
				int zeros=bonds-array.length;
				array=addZerosF(array,zeros);
				int[][] mat= new int[atomCount][atomCount];
				mat[0]=array;
				matrices.add(mat);
			}
		}
		List<int[][]> output= new ArrayList<int[][]>();
		generate(ac,1,6,matrices,output);
		if(verbose) {
			System.out.println("Number of hydrocarbons: "+output.size());
		}
		long endTime = System.nanoTime()- startTime;
        double seconds = (double) endTime / 1000000000.0;
		DecimalFormat d = new DecimalFormat(".###");
		if(verbose) {
			System.out.println("Duration:"+" "+d.format(seconds));
		}
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
				addBond(ac2,i,j,adjacency[i][j]);
			}
		}
		return ac2;
	}
	
	/**
	 * Add a bond to a atomcontainer
	 * @param ac atomcontainer
	 * @param i first index
	 * @param j second index
	 * @param order bond order
	 */
	
	public static void addBond(IAtomContainer ac, int i, int j, int order) {
		if(order==1) {
			ac.addBond(i, j, Order.SINGLE);
		}else if(order==2) {
			ac.addBond(i, j, Order.DOUBLE);
		}else if(order==3) {
			ac.addBond(i, j, Order.TRIPLE);
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
	
	public static void main(String[] arguments) throws FileNotFoundException, UnsupportedEncodingException, CloneNotSupportedException {
		HydrocarbonGenerator gen= new HydrocarbonGenerator();
		//String[] arguments1= {"-f","C6H12","-v"};
		try {
			gen.parseArguments(arguments);
			HydrocarbonGenerator.run(HydrocarbonGenerator.formula);
		} catch (Exception e) {
			if (HydrocarbonGenerator.verbose) e.getCause(); 
		}
	}
}
