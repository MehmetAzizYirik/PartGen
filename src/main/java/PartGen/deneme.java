package PartGen;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
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
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.graph.invariant.Canon;
import org.openscience.cdk.graph.invariant.MorganNumbersTools;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IIsotope;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.signature.MoleculeSignature;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;


public class deneme {
	public static IChemObjectBuilder builder =SilentChemObjectBuilder.getInstance();
	public static LinkedHashSet<int[][]> outputCheck = new LinkedHashSet<int[][]>();
	public static LinkedHashSet<int[]> arrayCheck = new LinkedHashSet<int[]>();
	public static IMolecularFormula formula=null;
	public static IAtomContainer acontainer;
	public static Map<Integer, Integer> capacities;
	public static boolean verbose = false;
	public static int size;
	public static int isotopes;
	public static int totalHydrogen;
	public static int[] capacity;
	public static int[] valences;
	public static int[] isotopesCount; // Total number of atoms.
	
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
	 * Add an element to an int array.
	 * @param a int array
	 * @param e int value to add
	 * @return
	 */
	
	public static int[] addElement(int[] a, int e) {
        a  = Arrays.copyOf(a, a.length + 1);
        a[a.length - 1] = e;
        return a;
    }
	
	/**
	 * Factorial of an integer.
	 * @param i integer
	 * @return
	 */
	
	public static int factorial(int i){
		if (i==0){
			return 1;
		}
		return i * factorial(i - 1);
	}
	
	/**
	 * Combination of two integers
	 * @param m integer 
	 * @param n integer
	 * @return
	 */
	
	public static int combination(int m, int n) {
		return factorial(m) / (factorial(n) * factorial(m - n));
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
	 * Sum all array entries.
	 * @param array int array
	 * @return
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
	 * For aligning two arrays. 
	 * @param arrays list of int arrays
	 * @return
	 */
	
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
	
	public static int[] arrayAlign(int[] a, int[] b) {
		List<int[]> arrays= new ArrayList<int[]>();
		arrays.add(a);
		arrays.add(b);
		return mergeArrays(arrays);
	}
	
	/**
	 * For the partitioning function, combining all the generated arrays. 
	 * @param lists list of int[] lists.
	 * @return
	 */
	
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
	            	newComb.add(arrayAlign(arr1,arr2));
	            }
	        }
	        comb = newComb;
	    }
	    return comb;
	}
	
	/**
	 * To set entries of a matrix. Both the row and their transpose columns. 
	 * @param mat int matrix
	 * @param index row index
	 * @param array entry values
	 * @return int matrix with the new values
	 */
	
	public static int[][] setEntries(int[][] mat,int index, int[] array) {
		for(int i=0;i<array.length;i++) {
			mat[((i+index)+1)][index]=array[i];
			mat[index][((i+index)+1)]=array[i];
		}
		return mat;
	}
	
	/**
	 * Setting all general values such as valences, isotopes and capacities.
	 * @param formula Input molecular formula. 
	 * @return
	 */
	
	public static int[] setValues(IMolecularFormula formula) {
		int i=0;
		int[] isotopes = new int[formula.getIsotopeCount()];
		for(IIsotope top:formula.isotopes()) {
			isotopes[i]=formula.getIsotopeCount(top);
			i++;
		}
		deneme.valences=valences(acontainer);
		deneme.isotopesCount=isotopes;
		return valences;
	}
	
	/**
	 * Returning a valence array for all the atoms in the molecule. 
	 * Useful for dealing with valences in row by row generation.
	 * 
	 * @param ac acontainer 
	 * @return
	 */
	
	public static int[] valences(IAtomContainer ac) {
		int[] valence = new int[size];
		int i=0;
		for(IAtom atom: ac.atoms()) {
			valence[i]=capacities.get(atom.getAtomicNumber())+1;
			i++;
		}
		return valence;
	}
	
	/**
	 * To find the starting index of the last isotope in an atom container.
	 * @return index 
	 */
	
	public static int lastIsotopeIndex() {
		String symbol=acontainer.getAtom(0).getSymbol();
		int index=0;
		for(int i=1;i<size;i++) {
			String symbol2=acontainer.getAtom(i).getSymbol();
			if(!symbol2.equals(symbol)) {
				symbol=symbol2;
				index=i;	
			}
		}
		return index;
	}
	
	/**
	 * Counting the occurrence of all isotopes in the row.
	 * @param index row index
	 * @return
	 */
	
	public static int[] isotopes(int index) {
		int[] isotope= new int[isotopes];
		String symbol=acontainer.getAtom(index+1).getSymbol();
		int j=0;
		isotope[j]++;
		for(int i=index+2;i<size;i++) {
			String symbol2=acontainer.getAtom(i).getSymbol();
			if(!symbol2.equals(symbol)) {
				symbol=symbol2;
				j++;
				isotope[j]++;
			}else {
				isotope[j]++;
			}
		}
		return isotope;
	}
	
	/**
	 * Counting number of isotopes for a given row index
	 * @param index row index in the matrix
	 * @return
	 */
	
	public static int countIsotopes(int index) {
		int[] iso = isotopes(index);
		int count=0;
		for(int i=0;i<iso.length;i++) {
			if(iso[i]!=0) {
				count++;
			}
		}
		return count;
	}

	/**
	 * When there are multiple isotopes, we use the limit. Untill the limit it should
	 * generate with multigen, then just simpleBondAdd for the remaining. 
	 * @param bonds number of bonds
	 * @param mat input int matrix
	 * @param limit the index where the last isotope type starts in the matrix
	 * @return
	 * @throws CloneNotSupportedException
	 */
	
	public static int[][] distribute1s(int bonds,int[][] mat,int limit) throws CloneNotSupportedException {
		int count=0;
		int[][] copy=copy(mat);
		for(int i=limit;i<size;i++) {
			for(int j=i+1;j<size;j++) {
				for(int k=0;k<=(valences[i]-sum(mat[i]));k++) {
					count++;
					if(count<=bonds) {
						copy[i][j]=1;
						copy[j][i]=1;
					}
				}
			}
		}
		return copy;
	}
	
	/**
	 * When there is only one isotope type, then we simply add bonds without any order. Since all 
	 * atoms are the same. 
	 * @param bonds
	 * @param mat
	 * @return
	 * @throws CloneNotSupportedException
	 */
	
	public static int[][] distribute1s(int bonds) throws CloneNotSupportedException {
		int count=0;
		int valenceValue= valences[0];
		int[][] mat= new int[size][size];
		for(int i=0;i<size;i++) {
			for(int j=i+1;j<size;j++) {
				for(int k=0;k<=(valenceValue-sum(mat[i]));k++) {
					count++;
					if(count<=bonds) {
						mat[i][j]=1;
						mat[j][i]=1;
					}
				}
			}
		}
		return mat;
	}
	
	/**
	 * Partitioning given number of bonds into arrays, might starts with 1 or 0. No strict rule.
	 * @param total number of bonds to distribute
	 * @param n		number of bonds to distribute
	 * @param d		the array size
	 * @param depth always starts with 0.
	 * @return
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
				if(i<=valences[item.length]) {
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
	
	/**
	 * Partitioning given number of bonds into arrays, always starting with 1.
	 * @param total number of bonds to distribute
	 * @param n		number of bonds to distribute
	 * @param d		the array size
	 * @param depth always starts with 0.
	 * @return
	 */
	
	public static List<int[]> partition1s(int total,int n, int d,int depth) {
		if(d==depth) {
			List<int[]> array= new ArrayList<int[]>();
			int[] take=new int[0];
			array.add(take);
			return array;
		}
		return buildArray1s(total,n,d,depth);
	}
	
	public static List<int[]> buildArray1s(int total,int n,int d, int depth){
		List<int[]> array= new ArrayList<int[]>();
		IntStream range = IntStream.rangeClosed(0,1);
		for(int i:range.toArray()) {
			for(int[] item: partition1s(total,n-i,d,depth+1)) {
				if(i<valences[item.length]) {
					item=addElement(item,i);
			        if(item.length==d) {
			        	if(sum(item)<=total) {
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
	
	/**
	 * For therow index, generating all possible bond distributions. These distributions are arrays
	 * filled with 1s and 0s.
	 * @param index row index
	 * @param mat adjacency matrix
	 * @return
	 * @throws CloneNotSupportedException
	 */
	
	public static List<int[]> rowsGen(int index, int[][] mat) throws CloneNotSupportedException{
		List<int[]> rows= new ArrayList<int[]>();
		if(countIsotopes(index)==1) {
			int value=Math.min(valences[index]-sum(mat[index]), (size-(index+1)));
			for(int[] array:partition1s(value,value,size-(index+1),0)) {
				rows.add(array);
			}
		}else {
			int[] iso= isotopes(index);	
			for(int[] arr:partition(valences[index],valences[index],countIsotopes(index),0)){
				LinkedList<List <int[]>> lists = new LinkedList<List <int[]>>();
				for(int i=0;i<arr.length;i++) {
					List<int[]> list=partition1s(arr[i],arr[i],iso[i],0); 
					lists.add(list);
				}
				List<int[]> combined=combineArrays(lists);
				for(int[] c:combined) {
					if (arrayCheck.add(c)) rows.add(c);
				}
			}
		}
		return rows;
	}
	
	/**
	 * For the formulas with more than 1 isotope types, first we need to generate all possible 
	 * 1 and 0 distributions in rows. Once the limit is reached, the simpleBondAdd is terminated.
	 * 
	 * @param index
	 * @param limit
	 * @param matrices
	 * @param output
	 * @return
	 * @throws CloneNotSupportedException
	 */
	
	public static Set<int[][]> multiGen(int index,int limit,Set<int[][]> matrices,Set<int[][]> output) throws CloneNotSupportedException {
		if(index==limit) {
			for(int[][] mat:matrices) {
				if(valenceCheck(mat)) { 
					if (outputCheck.add(mat)) output.add(mat);
				}
			}
		}else {
			if(index==0) {
				int[][] initial= new int[size][size];
				List<int[]> firstRows = rowsGen(0,initial); // to remove the duplicates I should add count down.
				for(int[] arr:firstRows){
					if(sum(arr)!=0) {
						int[][] mat= new int[size][size];
						setEntries(mat,index,arr);
						matrices.add(mat);
						multiGen((index+1),limit,matrices,output);
					}
				}
			}else {
				for(int[][] mat:matrices) {
					Set<int[][]> list= new HashSet<int[][]>();
					List<int[]> rows= rowsGen(index,mat);
					for(int[] arr:rows) {
						int[][] copy=copy(mat);
						setEntries(copy,index,arr);
						list.add(copy);
						multiGen((index+1),limit,list,output);
					}
				}
			}
		}
		return matrices;
	}
	
	/**
	 * Starting from the row, counting all the remaining open sites of atoms. 
	 * @param row row index
	 * @param mat int matrix
	 * @return
	 */
	
	public static int remaining(int row, int[][] mat) {
		int count=0; 
		for(int i=row;i<size;i++) {
			count=count+(valences[row]-sum(mat[row]));
		}
		return count;
	}
	
	/**
	 * Starting from the limit, simply adding the given number of bonds to the matrix.
	 * @param mat int matrix filled by multiGen function.
	 * @param limit the index where the last isotope starts in the matrix
	 * @return
	 * @throws CloneNotSupportedException
	 */
	
	public static List<int[][]> simpleBondAdd(int[][] mat, int limit) throws CloneNotSupportedException {
		int remaining=remaining(limit,mat);
		List<int[][]> matrices= new ArrayList<int[][]>();
		for(int i= remaining;i>=0;i--) { 
			int[][] mat2= distribute1s(i,mat,limit);
			if(zeroColumnCheckAll(mat2)) {
				if(sum(mat2)>=size-1 && valenceCheck(mat2)) {
					matrices.add(mat2); 
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
	 * ZeroColumnCheck in a matrix. If there is a zero column, then the corresponding atom does not have any interactions.
	 * @param mat int matrix
	 * @return
	 */
	
	public static boolean zeroColumnCheck(int[][] mat) {
		boolean check=true;
		if(sum(getColumn(mat,size-1))==0) {
			check=false;
		}
		return check;
	}
	//TODO: I dont know whether there should be a zero column or not.Should check.
	
	public static boolean zeroColumnCheckAll(int[][] mat) {
		boolean check=true;
		for(int i=0;i<size;i++) {
			if(sum(getColumn(mat,i))==0) {
				check=false;
			}
		}
		return check;
	}
	
	/**
	 * Generating atomContainers for a list of adjacency matrices.
	 * @param adcacencyMatrices int matrices
	 * @return
	 * @throws CloneNotSupportedException
	 * @throws CDKException 
	 */
	
	public static List<IAtomContainer> generateAtomContainers(Set<int[][]> adjacencyMatrices) throws CloneNotSupportedException, CDKException{
		Set<String> check = new HashSet<String>();
		ArrayList<IAtomContainer> acontainers = new ArrayList<IAtomContainer>();
		for(int[][] adjacency:adjacencyMatrices) {
			IAtomContainer ac=acontainer.clone();
			IAtomContainer newAc=setBonds(ac,adjacency);
			MoleculeSignature molSig = new MoleculeSignature(newAc);			
			String can=molSig.toCanonicalString();
			if (check.add(can)) {
				acontainers.add(newAc);
			}
		}
		return acontainers;
	}
	
	/**
	 * Setting bonds to an atomContainer
	 * @param ac atomContainer
	 * @param adjacency adjacency matrix
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
	
	/**
	 * Molecule Depiction
	 * @param mol atomContainer
	 * @param path part to generate png file
	 * @throws CloneNotSupportedException
	 * @throws CDKException
	 * @throws IOException
	 */
	
	public static void depict(IAtomContainer mol, String path) throws CloneNotSupportedException, CDKException, IOException{
		DepictionGenerator depict = new DepictionGenerator();
		depict.withCarbonSymbols().withSize(1000, 1000).withZoom(4).depict(mol).writeTo(path);
	}
			
	/**
	 * To initialise the inputs and run the functions while recording the duration time.
	 * @throws CDKException 
	 * @throws IOException 
	 */
	
	public static List<IAtomContainer> run(IMolecularFormula formula) throws CloneNotSupportedException, CDKException, IOException {
		long startTime = System.nanoTime(); //Recording the duration time.
		int hydrogen=formula.getIsotopeCount(builder.newInstance(IIsotope.class, "H"));
		String formulaString =MolecularFormulaManipulator.getString(formula);
		deneme.isotopes=formula.getIsotopeCount()-1;//If formula has H.
		formula.removeIsotope(builder.newInstance(IIsotope.class, "H"));
		IAtomContainer ac=MolecularFormulaManipulator.getAtomContainer(formula);
		deneme.size=ac.getAtomCount();
		deneme.acontainer=ac;
		deneme.totalHydrogen=hydrogen;
		setValues(formula);
		Set<int[][]> matrices = new HashSet<int[][]>();
		Set<int[][]> output = new HashSet<int[][]>();
		Set<int[][]> output2 = new HashSet<int[][]>();
		Set<int[][]> result= new HashSet<int[][]>();
		if(verbose) {
			System.out.println("For molecular formula "+ formulaString +", calculating all the possible distributions of "+totalHydrogen+" "+"hydrogens ..." );
		}
		int count=0;
		if(isotopes==1) {
			int valenceTotal= valences[0]*isotopesCount[0];
			for(int b= valenceTotal;b>size-1;b--) {
				output.add(distribute1s(b));
			}
			result= output;
		}else {
			int limit= lastIsotopeIndex(); 
			multiGen(0,limit,matrices,output);
			for(int[][] mat: output) {
				output2.addAll(simpleBondAdd(mat,limit));
			}
			result=output2;
		}
		long endTime = System.nanoTime()- startTime;
        double seconds = (double) endTime / 1000000000.0;
		DecimalFormat d = new DecimalFormat(".###");
		List<IAtomContainer> acontainers=generateAtomContainers(result);
		count=acontainers.size();
		if(verbose) {
			System.out.println("Number of distributions: "+count);
			System.out.println("Duration:"+" "+d.format(seconds));
		}
		return acontainers;
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
	public static void main(String[] arguments) throws FileNotFoundException, UnsupportedEncodingException, CloneNotSupportedException {
		deneme distribution= new deneme();
		String[] arguments1= {"-f","CCH4O","-v"};
		try {
			distribution.parseArguments(arguments1);
			int count=0;
			for(IAtomContainer ac:deneme.run(deneme.formula)) {
				deneme.depict(ac, "C:\\Users\\mehme\\Desktop\\output\\"+count+".png");
				count++;
			}
		} catch (Exception e) {
			if (deneme.verbose) e.getCause(); 
		}
	}
}
