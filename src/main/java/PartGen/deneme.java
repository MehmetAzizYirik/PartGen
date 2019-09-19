package PartGen;

import java.io.FileNotFoundException;
import java.io.UnsupportedEncodingException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
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


public class deneme {
	public static IChemObjectBuilder builder =SilentChemObjectBuilder.getInstance();
	public static IMolecularFormula formula=null;
	public static IAtomContainer acontainer;
	public static Map<Integer, Integer> capacities;
	public static boolean verbose = false;
	public static int size;
	public static int isotopes;
	public static int[] capacity;
	public static int[] valences;
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
		int[] capacity = new int[formula.getIsotopeCount()];
		int[] valences = new int[formula.getIsotopeCount()];
		int[] totalAtom = new int[formula.getIsotopeCount()];
		int i=0;
		for(IIsotope top:formula.isotopes()) {
			totalAtom[i]=formula.getIsotopeCount(top);
			valences[i]=capacities.get(top.getAtomicNumber());
			capacity[i]=capacities.get(top.getAtomicNumber())*formula.getIsotopeCount(top);
			i++;
		}
		deneme.capacity=capacity;
		deneme.valences=valences;
		deneme.totalAtom=totalAtom;
		return capacity;
	}
	
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
	public static List<int[]> combineArrays(LinkedList<List <int[]>> lists) {
		List<int[]> comb = new ArrayList<int[]>();
	    for (int[] s: lists.removeFirst()) {
	    	comb.add(s);
	    }
	    while (!lists.isEmpty()) {
	        List<int[]> list = lists.removeFirst();
	        List<int[]> newComb =  new ArrayList<int[]>();
	        for (int[] arr1: comb) { 
	            for (int[] arr2 : list) { 
	            	newComb.add(arraySum(arr1,arr2));
	            }
	        }
	        comb = newComb;
	    }
	    return comb;
	}
	
	public static List<int[]> alignArrays(LinkedList<List <int[]>> lists) {
		List<int[]> comb = new ArrayList<int[]>();
	    for (int[] s: lists.removeFirst()) {
	    	comb.add(s);
	    }
	    while (!lists.isEmpty()) {
	        List<int[]> list = lists.removeFirst();
	        List<int[]> newComb =  new ArrayList<int[]>();
	        for (int[] arr1: comb) { 
	            for (int[] arr2 : list) { 
	            	newComb.add(arraySum(arr1,arr2));
	            }
	        }
	        comb = newComb;
	    }
	    return comb;
	}
	/**
	 * To initialise the inputs and run the functions while recording the duration time.
	 * @throws CDKException 
	 */
	
	public static List<IAtomContainer> run(IMolecularFormula formula) throws FileNotFoundException, UnsupportedEncodingException, CloneNotSupportedException, CDKException {
		long startTime = System.nanoTime(); //Recording the duration time.
		int hydrogen=formula.getIsotopeCount(builder.newInstance(IIsotope.class, "H"));
		String formulaString =MolecularFormulaManipulator.getString(formula);
		deneme.isotopes=formula.getIsotopeCount()-1;
		formula.removeIsotope(builder.newInstance(IIsotope.class, "H"));
		IAtomContainer ac=MolecularFormulaManipulator.getAtomContainer(formula);
		deneme.size=ac.getAtomCount();
		deneme.acontainer=ac;
		setValues(formula);
		if(verbose) {
			System.out.println("For molecular formula "+ formulaString +", calculating all the possible distributions of "+totalHydrogen+" "+"hydrogens ..." );
		}
		List<int[]> result= new ArrayList<int[]>();
		int count=0;
		if(isotopes==1) {
			List<int[]> iarrays= new ArrayList<int[]>();
			int[] array = new int[0];
			distribute(iarrays,array,valences[0],totalAtom[0]);
			count=iarrays.size();
			result= iarrays;
		}else {
			List<int[]> distributions= new ArrayList<int[]>();
			for(int[] dene:partition(totalHydrogen,isotopes,0)){
				LinkedList<List <int[]>> lists = new LinkedList<List <int[]>>();
				for(int i=0;i<dene.length;i++) {
					deneme.hydrogens2distribute=dene[i];
					List<int[]> iarrays= new ArrayList<int[]>();
					int[] array = new int[0];
					distribute(iarrays,dene[i],array,valences[i],totalAtom[i]);
					lists.add(iarrays);
				}	
				List<int[]> combined=combineArrays(lists);
				count+=combined.size();
				distributions.addAll(combined);
			}
			result=distributions;
		}
		if(verbose) {
			System.out.println("Number of distributions: "+count);
		}
		long endTime = System.nanoTime()- startTime;
        double seconds = (double) endTime / 1000000000.0;
		DecimalFormat d = new DecimalFormat(".###");
		if(verbose) {
			System.out.println("Duration:"+" "+d.format(seconds));
		}

		return generateAtomContainers(result);
	}
	
	/**
	 * These functions are built for the integer partitioning problem.
	 */
	
	public static List<int[]> partition(int total,int n, int d,int depth) {
		if(d==depth) {
			List<int[]> array= new ArrayList<int[]>();
			int[] take=new int[0];
			array.add(take);
			return array;
		}
		return buildArray(total,n,d,depth);
		
	}
	
	public static List<int[]> buildArray(int total,int n,int d, int depth){
		List<int[]> array= new ArrayList<int[]>();
		IntStream range = IntStream.rangeClosed(0,n);
		for(int i:range.toArray()) {
			for(int[] item: partition(total,n-i,d,depth+1)) {
				if(i<=capacity[item.length]) {
					item=addElement(item,i);
					if(item.length==d) {
						if(sum(item)==total) {
							array.add(item);
						}
					}else {
						array.add(item);
					}
				}
			}
		}
		return array;
	}
	
	public static List<int[][]> complexGen(int index,List<int[][]> matrices,List<int[][]> output) throws CloneNotSupportedException {
		if(index==size-1) {
			for(int[][] mat:matrices) {
				if(valenceCheck(mat) && sum(mat)>=size-1) { //Minimal number of bonds is the number of atoms minus 1.
					//System.out.println(Arrays.deepToString(mat));
					output.add(mat);
				}
			}
		}else {
			for(int[][] mat:matrices) {
				List<int[]> iarrays= new ArrayList<int[]>();
				List<int[][]> list= new ArrayList<int[][]>();
				List<int[]> rows= new ArrayList<int[]>();
				for(int[] arr:partition(size-(index+1),size-(index+1),size-(index+1),0)){
					LinkedList<List <int[]>> lists = new LinkedList<List <int[]>>();
					for(int i=0;i<arr.length;i++) {
						lists.add(distribute1s(arr[i],capacity[i]));
					}
					List<int[]> combined=combineArrays(lists);
					rows.addAll(combined);
				}
				for(int[] arr:rows) {
					int zeros=size-arr.length;
					arr=addZerosF(arr,zeros);
					int[][] copy=copy(mat);
					copy[index]=arr;
					list.add(copy);
					complexGen(index+1,list,output);
				}
			}
		}
		return matrices;
	}
	
	public static List<int[][]> simpleGen(int index,List<int[][]> matrices,List<int[][]> output) throws CloneNotSupportedException{
		if(index==size-1) {
			for(int[][] mat:matrices) {
				if(valenceCheck(mat) && sum(mat)>=size-1) { 
					output.add(mat);
				}
			}
		}else {
			for(int[][] mat:matrices) {
				List<int[][]> list= new ArrayList<int[][]>();
				for(int[] arr: distribute1s(size-(index+1),size-(index+1))) {
					int zeros=size-arr.length;
					arr=addZerosF(arr,zeros);
					int[][] copy=copy(mat);
					copy[index]=arr;
					list.add(copy);
					simpleGen(index+1,list,output);
				}
			}
		}
		return matrices;
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
	
	public static int[] addZeros(int[] array, int zeros) {
		for(int i=0;i<zeros;i++) {
			array=addElement(array,0);
		}
		return array;
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
	
	public static List<int[]> distribute(List<int[]> arrays,int[]arr,int valence, int numAtom) throws CloneNotSupportedException {
		if(hydrogen==0 && sum(arr)==hydrogens2distribute){
			if(arr.length!=numAtom) {
				arr=addZeros(arr,(numAtom-arr.length));
			}
			arrays.add(arr);
		}else if((numAtom-arr.length)==1) {
			int add=Math.min(hydrogen,valence);
			hydrogen=hydrogen-add;
			if(arr.length==0) {
				distribute(arrays,0,addElement(arr,add),valence,numAtom); 
			}
			if((arr.length)>0) {
				if(arr[arr.length-1]<=add){ 
					distribute(arrays,0,addElement(arr,add),valence,numAtom);
				}
			}
		}else { 
			for(int i = Math.min(valence,hydrogen); i > 0; i--) {
				if(arr.length==0) {
					distribute(arrays,(hydrogen-i),addElement(arr,i),valence,numAtom); 
				}
				if((arr.length)>0) {
					if(arr[arr.length-1]<=i){ 
						distribute(arrays,(hydrogen-i),addElement(arr,i),valence,numAtom);
					}
				}
			}
		}
		return arrays;
	}
	
	public static List<int[]> distribute1s(int dis, int space) throws CloneNotSupportedException {
		List<int[]> arrays= new ArrayList<int[]>();
		for(int i = dis; i > 0; i--) {
			int[] arr= new int[space];
			arrays.add(add1s(arr,i));
		}
		return arrays;
	}
	
	public static int[] add1s(int[] array, int count) {
		int size= array.length;
		for(int i=size;i<count+size;i++) {
			array[i]=1;
		}
		return array;
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
			deneme.formula=MolecularFormulaManipulator.getMolecularFormula(formula, builder);
			if (cmd.hasOption("verbose")) deneme.verbose = true;
		
		} catch (ParseException e) {
			HelpFormatter formatter = new HelpFormatter();
			formatter.setOptionComparator(null);
			String header = "\nFor a molecular formula, it calculates all the possible hydrogen distributions to the atoms.";
			String footer = "\nPlease report issues at https://github.com/MehmetAzizYirik/deneme";
			formatter.printHelp( "java -jar deneme.jar", header, options, footer, true );
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
	public static void main(String[] arguments) throws FileNotFoundException, UnsupportedEncodingException {
		
		/**deneme distribution= new deneme();
		//String[] arguments1= {"-f","C7H8ClNO2S","-v"};
		try {
			distribution.parseArguments(arguments);
			deneme.run(deneme.formula);
		} catch (Exception e) {
			if (deneme.verbose) e.getCause(); 
		}**/
	}
}
