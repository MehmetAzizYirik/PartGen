package PartGen;

import java.io.FileNotFoundException;
import java.io.IOException;
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
import org.openscience.cdk.depict.DepictionGenerator;
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


public class deneme {
	public static IChemObjectBuilder builder =SilentChemObjectBuilder.getInstance();
	public static IMolecularFormula formula=null;
	public static IAtomContainer acontainer;
	public static Map<Integer, Integer> capacities;
	public static boolean verbose = false;
	public static int size;
	public static int isotopes;
	public static int totalHydrogen;
	public static int[] capacity;
	public static int[] valence;
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
		deneme.valence=valence(acontainer);
		deneme.totalAtom=totalAtom;
		return capacity;
	}
	
	public static int[] valence(IAtomContainer ac) {
		int[] valence = new int[size];
		int i=0;
		for(IAtom atom: ac.atoms()) {
			valence[i]=capacities.get(atom.getAtomicNumber())+1;
			i++;
		}
		return valence;
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
	
	public static int newIsotope() {
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
	
	public static boolean zeroColumnCheck(int[][] mat) {
		boolean check=true;
		if(sum(getColumn(mat,size-1))==0) {
			check=false;
		}
		return check;
	}
	/**
	 * To initialise the inputs and run the functions while recording the duration time.
	 * @throws CDKException 
	 * @throws IOException 
	 */
	
	public static void run(IMolecularFormula formula) throws CloneNotSupportedException, CDKException, IOException {
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
		int limit= newIsotope();
		int remaining=remaining(limit);
		int say=2;
		for(int[] r:rowsGen2(0)) {
			if(sum(r)!=0) {
				if(sum(r)==say) {
					System.out.println(Arrays.toString(r)+" "+say);
					int zeros=size-r.length;
					r=addZerosF(r,zeros);
					int[][] mat= new int[size][size];
					mat[0]=r;
					int[][] mat2=distribute1s(remaining,mat,limit);
					if(zeroColumnCheck(mat2)) {
						output2.add(mat2);
					}
					say--;
				}
			}
		}

		List<int[]> firstRows = rowsGen(0);
		for(int[] arr:firstRows){
			if(sum(arr)!=0) {
				int zeros=size-arr.length;
				arr=addZerosF(arr,zeros);
				int[][] mat= new int[size][size];
				mat[0]=arr;
				matrices.add(mat);
			}
		}
		
		multiGen(1,limit,matrices,output);
		for(int[][] mat: output) {
			output2.addAll(simpleBondAdd(ac,mat,limit));
		}
		//count=output.size();
		//result=output;
		int count2=0;
		for(int[][] mat: output2) {
			System.out.println(count2+" "+Arrays.deepToString(mat));
			count2++;
		}
		int count=0;
		for(IAtomContainer acon: generateAtomContainers(output2)) {
			depict(acon,"C:\\Users\\mehme\\Desktop\\output\\"+count+".png");
			count++;
		}
		
		/**
		setValues(formula);
		if(verbose) {
			System.out.println("For molecular formula "+ formulaString +", calculating all the possible distributions of "+totalHydrogen+" "+"hydrogens ..." );
		}
		Set<int[][]> result= new HashSet<int[][]>();
		int count=0;
		if(isotopes==1) {
			Set<int[][]> matrices= new HashSet<int[][]>();
			Set<int[][]> output= new HashSet<int[][]>();
			for(int[] array:distribute1s(size-1,size-1)) {
				if(sum(array)!=0) {
					int zeros=size-array.length;
					array=addZerosF(array,zeros);
					int[][] mat= new int[size][size];
					mat[0]=array;
					matrices.add(mat);
				}
			}
			
			simpleGen(1,matrices,output);
			count=output.size();
			result= output;
			int i=0;
			for(int[][] r: result) {
				System.out.println(Arrays.deepToString(r));
			}
			for(IAtomContainer acon: generateAtomContainers(result)) {
				depict(acon,"C:\\Users\\mehme\\Desktop\\output\\"+i+".png");
				i++;
			}
		}else {
			Set<int[][]> matrices = new HashSet<int[][]>();
			Set<int[][]> output = new HashSet<int[][]>();
			List<int[]> firstRows = rowsGen(0);
			for(int[] arr:firstRows){
				if(sum(arr)!=0) {
					int zeros=size-arr.length;
					arr=addZerosF(arr,zeros);
					int[][] mat= new int[size][size];
					mat[0]=arr;
					matrices.add(mat);
				}
			}
			multiGen(1,matrices,output);
			count=output.size();
			result=output;
			for(int[][] mat: output) {
				System.out.println(Arrays.deepToString(mat));
			}
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
		int i=0;
		for(IAtomContainer acon: generateAtomContainers(result)) {
			depict(acon,"C:\\Users\\mehme\\Desktop\\output\\"+i+".png");
			i++;
		}
		return result;**/
	}
	
	// Molecule depiction Generator
	public static void depict(IAtomContainer mol, String path) throws CloneNotSupportedException, CDKException, IOException{
		DepictionGenerator depict = new DepictionGenerator();
		depict.withCarbonSymbols().withSize(1000, 1000).withZoom(4).depict(mol).writeTo(path);
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
				if(i<=valence[item.length]) {
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
	
	public static Set<int[][]> multiGen(int index,Set<int[][]> matrices,Set<int[][]> output) throws CloneNotSupportedException {
		if(index==size-1) {
			for(int[][] mat:matrices) {
				if(valenceCheck(mat) && sum(mat)>=size-1) { //Minimal number of bonds is the number of atoms minus 1.
					output.add(mat);
				}
			}
		}else {
			for(int[][] mat:matrices) {
				Set<int[][]> list= new HashSet<int[][]>();
				List<int[]> rows= rowsGen(index);
				for(int[] arr:rows) {
					int zeros=size-arr.length;
					arr=addZerosF(arr,zeros);
					int[][] copy=copy(mat);
					copy[index]=arr;
					list.add(copy);
					multiGen(index+1,list,output);
				}
			}
		}
		return matrices;
	}
	
	public static Set<int[][]> multiGen(int index,int limit,Set<int[][]> matrices,Set<int[][]> output) throws CloneNotSupportedException {
		if(index==limit) {
			for(int[][] mat:matrices) {
				if(valenceCheck(mat)) { //Minimal number of bonds is the number of atoms minus 1.
					output.add(mat);
				}
			}
		}else {
			for(int[][] mat:matrices) {
				Set<int[][]> list= new HashSet<int[][]>();
				List<int[]> rows= rowsGen2(index);
				for(int[] arr:rows) {
					int[][] copy=copy(mat);
					copy[index]=arr;
					list.add(copy);
					multiGen(index+1,limit,list,output);
				}
			}
		}
		return matrices;
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
	
	public static List<int[][]> simpleGenBond(IAtomContainer ac) throws CloneNotSupportedException {
		int size= ac.getAtomCount();
		int maxBonds=combination(size,2);
		List<int[][]> matrices= new ArrayList<int[][]>();
		for(int i= maxBonds;i>size-2;i--) {
			matrices.add(distribute1s(i));
		}
		int[][] array= new int[size][size];
		for(int i=1;i<size-1;i++) {
			array[0][i]=1;
		}
		array[2][size-1]=1;
		matrices.add(array);
		return matrices;
	}
	
	public static List<int[][]> simpleBondAdd(IAtomContainer ac, int[][] mat, int limit) throws CloneNotSupportedException {
		int size= ac.getAtomCount();
		int remaining=remaining(limit);
		List<int[][]> matrices= new ArrayList<int[][]>();
		for(int i= remaining-1;i>=(size-1-sum(mat));i--) {
			int[][] mat2= distribute1s(i,mat,limit);
			if(zeroColumnCheck(mat2)) {
				matrices.add(distribute1s(i,mat,limit));
			}
		}
		/**int[][] array= new int[size][size];
		for(int i=1;i<size-1;i++) {
			array[0][i]=1;
		}
		array[2][size-1]=1;
		if(sum(array)>=size-1 && valenceCheck(mat)) {
			matrices.add(array);
		}**/
		return matrices;
	}
	
	public static Set<int[][]> simpleGen(int index,Set<int[][]> matrices,Set<int[][]> output) throws CloneNotSupportedException{
		if(index==size-1) {
			for(int[][] mat:matrices) {
				if(valenceCheck(mat) && sum(mat)>=size-1) { 
					output.add(mat);
				}
			}
		}else {
			for(int[][] mat:matrices) {
				Set<int[][]> list= new HashSet<int[][]>();
				for(int[] arr: distribute1s(size-(index+1),size-(index+1))) {
					//if(sum(arr)!=0) {
						int zeros=size-arr.length;
						arr=addZerosF(arr,zeros);
						int[][] copy=copy(mat);
						copy[index]=arr;
						list.add(copy);
					//}
				}
				simpleGen(index+1,list,output);
			}
		}
		return matrices;
	}
	
	public static List<int[]> rowsGen(int index) throws CloneNotSupportedException{
		List<int[]> rows= new ArrayList<int[]>();
		if(countIsotopes(index)==1) {
			int value=Math.min(valence[index], size-index-1);
			for(int[] array:partition1s(value,value,size-(index+1),0)) {
				int zeros=size-array.length;
				array=addZerosF(array,zeros);
				rows.add(array);
			}
		}else {
			int[] iso= isotopes(index);			
			for(int[] arr:partition(valence[index],valence[index],countIsotopes(index),0)){
				LinkedList<List <int[]>> lists = new LinkedList<List <int[]>>();
				for(int i=0;i<arr.length;i++) {
					List<int[]> list=partition1s(arr[i],arr[i],iso[i],0); 
					lists.add(list);
				}
				List<int[]> combined=combineArrays(lists);
				rows.addAll(combined);
			}
		}
		return rows;
	}
	public static List<int[]> rowsGen2(int index) throws CloneNotSupportedException{
		List<int[]> rows= new ArrayList<int[]>();
		int value=Math.min(valence[index], size-index-1);
		for(int v=value;v>0;v--) {
			for(int[] array:part1s(v,v,size-(index+1),0)) {
				int zeros=size-array.length;
				array=addZerosF(array,zeros);
				rows.add(array);
			}
		}
		return rows;
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
			if(!(sum(mat[i])<=valence[i] && sum(getColumn(mat,i))<=valence[i])) {
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
		int[] zeros= new int[space];
		arrays.add(zeros);
		for(int i = dis; i > 0; i--) {
			int[] arr= new int[space];
			for(int j=0;j<i;j++) {
				arr[j]=1;
			}
			arrays.add(arr);
		}
		return arrays;
	}
	
	/**
	 * These functions are built for the integer partitioning problem.
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
	
	public static List<int[]> part1s(int total,int n, int d,int depth) {
		if(d==depth) {
			List<int[]> array= new ArrayList<int[]>();
			int[] take=new int[0];
			array.add(take);
			return array;
		}
		return buildArr1s(total,n,d,depth);
	}
	
	public static List<int[]> buildArr1s(int total,int n,int d, int depth){
		List<int[]> array= new ArrayList<int[]>();
		IntStream range = IntStream.rangeClosed(0,1);
		for(int i:range.toArray()) {
			for(int[] item: part1s(total,n-i,d,depth+1)) {
				if(item.length==0) {
					if(array.size()==0) {
						item=addElement(item,1);
						array.add(item);
					}else {
						continue;
					}
				}else {
					if(i<2) {
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
		}
		return array;
	}
	public static List<int[]> buildArray1s(int total,int n,int d, int depth){
		List<int[]> array= new ArrayList<int[]>();
		IntStream range = IntStream.rangeClosed(0,1);
		for(int i:range.toArray()) {
			for(int[] item: partition1s(total,n-i,d,depth+1)) {
				if(i<valence[item.length]) {
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
	
	public static List<int[]> dist1s(int number, int index){
		List<int[]> arrays= new ArrayList<int[]>();
		IntStream range = IntStream.rangeClosed(0,1);
		int[] array= new int[size];
		for(int i=0;i<index+1;i++) {
			array[i]=0;
		}
		for(int i=0;i<size-index-1;i++) {
			
		}
		return arrays;
	}
	public static int remaining(int row) {
		int count=0;
		for(int i=row;i<size-1;i++) {
			count=count+(size-i-1);
		}
		return count;
	}

	public static int[][] distribute1s(int bonds,int[][] mat,int limit) throws CloneNotSupportedException {
		int count=0;
		int[][] copy=copy(mat);
		for(int i=limit;i<size;i++) {
			for(int j=i+1;j<size;j++) {
				count++;
				if(count<=bonds) {
					copy[i][j]=1;
				}
			}
		}
		return copy;
	}
	
	public static int[] add1s(int[] array, int count) {
		int size= array.length;
		int[] output= new int[size+count];
		for(int i=0;i<size;i++) {
			output[i]=array[i];
		}
		for(int i=size;i<count+size;i++) {
			output[i]=1;
		}
		return output;
	}
	
	/**
	 * Functions seting the implicit hydrogens of the atoms.
	 * @throws CloneNotSupportedException 
	 */
	
	public static List<IAtomContainer> generateAtomContainers(Set<int[][]> adcacencyMatrices) throws CloneNotSupportedException{
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
	public static void main(String[] arguments) throws FileNotFoundException, UnsupportedEncodingException, CloneNotSupportedException {
		/**int[] cp= new int[3];
		List<int[][]> matrices= new ArrayList<int[][]>();
		List<int[][]> output= new ArrayList<int[][]>();
		for(int[] arr:partition(8,8,3,0)){
			System.out.println(Arrays.toString(arr));
		}
		for(int[] array:distribute1s(3,3)) {
			int zeros=4-array.length;
			array=addZerosF(array,zeros);
			//TODO: first line can be zero 
			int[][] mat= new int[4][4];
			mat[0]=array;
			matrices.add(mat);
		}
		simpleGen(1,matrices,output);
		for(int[][] mat: output) {
			System.out.println(Arrays.deepToString(mat));
		}**/
	
		deneme distribution= new deneme();
		String[] arguments1= {"-f","C3H3OO","-v"};
		try {
			distribution.parseArguments(arguments1);
			deneme.run(deneme.formula);
		} catch (Exception e) {
			if (deneme.verbose) e.getCause(); 
		}
	}
}
