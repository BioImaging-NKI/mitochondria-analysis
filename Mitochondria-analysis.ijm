#@ File[] (label = "Input image stack", style = "file") inputFiles
#@ File (label = "Output folder", style = "directory") outputFolder
#@ Double (label = "Normalization lower percentage", value = 25, style="format:0.00") lowerPercentile
#@ Double (label = "Normalization upper percentage", value = 98, style="format:0.00") upperPercentile
#@ String thresholdMethod (label = "Threshold method for segmentation", choices={"Default", "Huang", "Intermodes", "IsoData", "IJ_IsoData", "Li", "MaxEntropy", "Mean", "MinError", "Minimum", "Moments", "Otsu", "Percentile", "RenyiEntropy", "Shanbhag", "Triangle", "Yen"}, value = "Otsu", style="listBox")
#@ Integer	dilateLabelsRadius (value=4, min=0, label="Connect structures within distance (pixels)")
#@ Integer (label = "Overlay opacity (%)", value = 50) opacity

ignoreZeros = true;
dilateLabelsRadius = dilateLabelsRadius/2; //because dilation goes in two directions

for (f = 0; f < inputFiles.length; f++) {
	
//	Overlay.remove;
	run("Close All");
	setOption("ScaleConversions", true);
	run("CLIJ2 Macro Extensions", "cl_device=");
	Ext.CLIJ2_clear();
	run("Select None");
	print("\\Clear");
	roiManager("reset");

	open(inputFiles[f]);
	open(File.getDirectory(inputFiles[f]) + File.getNameWithoutExtension(inputFiles[f]) + "_ROIs.zip");

	original = getTitle();
	getDimensions(width, height, channels, slices, frames);
	getVoxelSize(pw, ph, pd, unit);

	resultsTable = "Results per cel";
	if(!isOpen(resultsTable)) Table.create(resultsTable);
	else Table.reset(resultsTable);
	
	//3D segmentation
	//run("Voronoi Threshold Labler (2D/3D)", "filtermethod=DoG filterradius=1.0 backgroundsubtractionmethod=None backgroundradius=0.0 histogramusage=full thresholdmethod=MaxEntropy fillholes=Off separationmethod=[Eroded sphere] spotsigma=1.0 maximaradius=0.0 volumerange=100-Infinity excludeonedges=false outputtype=Labels stackslice=17 applyoncompleteimage=false processonthefly=false");
	//run("Voronoi Threshold Labler (2D/3D)", "filtermethod=DoG filterradius=0.5 backgroundsubtractionmethod=None backgroundradius=0.0 histogramusage=full thresholdmethod=MaxEntropy fillholes=Off separationmethod=None spotsigma=1.0 maximaradius=0.0 volumerange=100-Infinity excludeonedges=false outputtype=Labels stackslice=17 applyoncompleteimage=false processonthefly=false");
	
	//selectWindow("VTL_"+original);
	Stack.setSlice(floor(slices/2));
	//setBatchMode(true);
	
	nrCells = roiManager("count");
	nrOfStructures = newArray(nrCells);
	
	//newImage("Merged_labels", "16-bit black", width, height, slices);
	//run("glasbey_on_dark");
	//setBatchMode("show");
	
	selectImage(original);
	getDimensions(width, height, channels, slices, frames);
	//Remove 'black' holes (deconvolution artifacts)
	for (z=1; z<=slices; z++) {
		Stack.setSlice(z);
		changeValues(0, 100, 100);
	}
	
	setBatchMode(true);
	for (i = 0; i < nrCells; i++) {
		showStatus("Analyzing cell "+i+1+" / "+nrCells);
		showProgress(i, nrCells);
		
		selectImage(original);
		roiManager("select", i);
		getSelectionBounds(x, y, selectionWidth, selectionHeight);
		run("Duplicate...", "title=cell_"+i+1+" duplicate");
		run("Clear Outside", "stack");
		Overlay.remove;
		Stack.setSlice(floor(slices/2));
		
		//Set a plausible Min and Max value, required because VTL converts to 8-bit with scaling.
		run("Z Project...", "projection=[Max Intensity]");
		run("Restore Selection");
		minAndMax = get_percentile_values(lowerPercentile/100, upperPercentile/100, ignoreZeros);
		close();
		selectImage("cell_"+i+1);
		setMinAndMax(minAndMax[0], minAndMax[1]);
		print("cell "+i+1+" min and max: "+minAndMax[0]+", "+minAndMax[1]);
		setBatchMode("exit and display");
		run("Voronoi Threshold Labler (2D/3D)", "filtermethod=Gaussian filterradius=0.5 backgroundsubtractionmethod=None backgroundradius=0.0 histogramusage=[ignore black] thresholdmethod="+thresholdMethod+" fillholes=Off separationmethod=None spotsigma=1.0 maximaradius=0.0 volumerange=100-Infinity excludeonedges=false outputtype=Labels stackslice=16 applyoncompleteimage=false processonthefly=false");
		setBatchMode(true);
		labelmap_cell = getTitle();
	//	Ext.CLIJ2_voronoiOtsuLabeling(Image_input, Image_destination, Number_spot_sigma, Number_outline_sigma);
	
		run("Duplicate...", "title=cell_"+i+1+"_skeleton duplicate");
		skeleton_cell = getTitle();
		setThreshold(1, 65535);
		run("Convert to Mask", "black");
		run("Red");
		selectImage(skeleton_cell);
		run("Skeletonize (2D/3D)");
	
		selectImage(skeleton_cell);
		run("Analyze Skeleton (2D/3D)", "prune=none");
		taggedSkeleton = getTitle();
		run("Glasbey on dark");
		setMinAndMax(0, 197);
	
		branches = Table.getColumn("# Branches", "Results");
		branchlength = Table.getColumn("Average Branch Length", "Results");
		junctions = Table.getColumn("# Junctions", "Results");
		triple = Table.getColumn("# Triple points", "Results");
		quadruple = Table.getColumn("# Quadruple points", "Results");
	
		branches_total = sumArray(branches);
		branchlength_total = sumArray(branchlength);
		junctions_total = sumArray(junctions);
		triple_total = sumArray(triple);
		quadruple_total = sumArray(quadruple);
		
		network_length = 0;
		for (k = 0; k < nResults; k++) network_length += (branches[k] * branchlength[k]);
		average_branch_length = network_length/branches_total;
		
		//make labelmap isotropic
		selectImage(labelmap_cell);
		run("Scale...", "x=1.0 y=1.0 z="+pd/pw+" interpolation=None process create");
		labelmap_cell_isotropic = getTitle();

		selectImage(labelmap_cell_isotropic);
		getDimensions(width, height, channels, slices_isotropic, frames);
		Stack.setSlice(floor(slices_isotropic/2));
		Ext.CLIJ2_push(labelmap_cell_isotropic);
		close(labelmap_cell_isotropic);
		Ext.CLIJ2_threshold(labelmap_cell_isotropic, labelmap_cell_binary, 1);
	
		//Ext.CLIJ2_drawMeshBetweenProximalLabels(labels_isotropic, labels_mesh, dilateLabelsRadius*10);
		//Ext.CLIJ2_drawMeshBetweenTouchingLabels(Image_input, Image_destination);
		
		Ext.CLIJ2_dilateLabels(labelmap_cell_isotropic, labelmap_cell_dilated, dilateLabelsRadius);
		Ext.CLIJ2_release(labelmap_cell_isotropic);
	
	//Ext.CLIJ2_pull(labelmap_cell_dilated);
	//setBatchMode("show");
	//waitForUser(2);
	//close();
		Ext.CLIJ2_mergeTouchingLabels(labelmap_cell_dilated, labelmap_cell_dilated_merged);
		Ext.CLIJ2_release(labelmap_cell_dilated);
	
		Ext.CLIJ2_multiplyImages(labelmap_cell_dilated_merged, labelmap_cell_binary, labelmap_cell_merged);
		Ext.CLIJ2_release(labelmap_cell_dilated_merged);
		Ext.CLIJ2_release(labelmap_cell_binary);
	
		Ext.CLIJ2_closeIndexGapsInLabelMap(labelmap_cell_merged, labelmap_cell_merged_gapsclosed);
		Ext.CLIJ2_release(labelmap_cell_merged);
		Ext.CLIJ2_pull(labelmap_cell_merged_gapsclosed);
		Ext.CLIJ2_release(labelmap_cell_merged_gapsclosed);
		
		selectImage(labelmap_cell_merged_gapsclosed);
	//	setBatchMode("show");
	
		//De-isotropify
		run("Scale...", "x=1.0 y=1.0 z="+pw/pd+" interpolation=None process create");
		labelmap_cell_merged = "cell_"+i+1+"_labelmap";
		rename(labelmap_cell_merged);
		run("glasbey_on_dark");
		setMinAndMax(0, 255);
		close(labelmap_cell_merged_gapsclosed);
	
		Stack.getStatistics(voxelCount, mean, min, max, stdDev);
		nrOfStructures[i] = max;
	
		//Create mask and measure volume and intensity
		selectImage(labelmap_cell);
		setThreshold(1, 65535);
		run("Convert to Mask", "black");
		run("Intensity Measurements 2D/3D", "input=cell_"+i+1+" labels="+labelmap_cell+" mean stddev max min median mode skewness kurtosis numberofvoxels volume");
		mean_intensity = Table.get("Mean", 0, "cell_"+i+1+"-intensity-measurements");
		volume = Table.get("Volume", 0, "cell_"+i+1+"-intensity-measurements");

		//Fill resultsTable
		Table.set("cell nr", i, i+1, resultsTable);
		Table.set("nr of structures", i, nrOfStructures[i], resultsTable);
		Table.set("total length  ("+unit+")", i, network_length, resultsTable);
		Table.set("volume ("+unit+"^3)", i, volume, resultsTable);
		Table.set("mean intensity", i, mean_intensity, resultsTable);
		Table.set("total intensity", i, mean_intensity * volume, resultsTable);
		Table.set("intensity per running "+unit, i, mean_intensity*volume/network_length, resultsTable);
		Table.set("mean diameter ("+unit+")", i, 2*sqrt((volume/network_length)/PI), resultsTable);
		Table.set("average branch length", i, average_branch_length, resultsTable);
		Table.set("junctions", i, junctions_total, resultsTable);
		Table.set("# branch points", i, branches_total, resultsTable);
		Table.set("# triple points", i, triple_total, resultsTable);
		Table.set("# quadruple points", i, quadruple_total, resultsTable);
		Table.update(resultsTable);
	
		//Fix shifted labels (by one slice)
		selectImage(labelmap_cell_merged);
		Stack.setSlice(slices);
		run("Add Slice");
		Stack.setSlice(1);
		run("Delete Slice");
	
		//Overlay segmentations and tagged skeleton
		overlay_image_3D(original, labelmap_cell_merged, x, y, 30);
		overlay_image_3D(original, taggedSkeleton, x, y, 100);
	
	//	//Overlay the labels onto the original image
	//	for (z=1; z<=slices; z++) {
	//		selectImage("cell_"+i+1+"_labelmap");
	//		Stack.setSlice(z);
	//		selectImage(original);
	//		Stack.setSlice(z);
	//		run("Add Image...", "image=[cell_"+i+1+"_labelmap] x="+x+" y="+y+" opacity="+opacity+" zero");
	//	}
	
		close("cell_"+i+1);
		close(labelmap_cell);
		close(skeleton_cell);
		close(taggedSkeleton);
		close(labelmap_cell_merged);
		close("cell_"+i+1+"-intensity-measurements");
	
	}
	
	cellNrArray = addScalarToArray(Array.getSequence(nrCells),1);

	selectImage(original);
	saveAs("tiff", outputFolder + File.separator + File.getNameWithoutExtension(inputFiles[f]));
	selectWindow(resultsTable);
	Table.save(outputFolder + File.separator + File.getNameWithoutExtension(inputFiles[f]) + "_results.tsv");
}	


//Overlays a 3D image
function overlay_image_3D(image, overlay, x, y, opacity) {
	batchMode = is("Batch Mode");
	if(!batchMode) setBatchMode(true);
	selectImage(image);
	getDimensions(width, height, channels, slices, frames);
	Stack.getPosition(channel, slice, frame);
//	Overlay.clear;
	for (i = 1; i <= slices; i++) {
		selectImage(overlay);
		Stack.setSlice(i);
		selectImage(image);
		Stack.setSlice(i);
		run("Add Image...", "image=["+overlay+"] x="+x+" y="+y+" opacity="+opacity+" zero");
	}
	Stack.setPosition(channel, slice, frame);
	if(!batchMode) setBatchMode(false);
}


//Returns the sum of all elements of an arrays, ignoring NaNs
function sumArray(array) {
	sum=0;
	for (a=0; a<lengthOf(array); a++) {
		if(!isNaN(array[a])) sum=sum+array[a];
	}
	return sum;
}


//Returns the mean of the array
function meanOfArray(array) {
	Array.getStatistics(array, min, max, mean, stdDev);
	return mean;
}


//Returns the mean of all elements of an array; NaNs are ignored
function meanOfArrayExcludingNaNs(array) {
	sum=0;
	nans=0;
	for (a=0; a<lengthOf(array); a++) {
		if(!isNaN(array[a])) sum=sum+array[a];
		else nans+=1;
	}
	return sum/(array.length-nans);
}


//Returns the stdDev of the array
function stdDevOfArray(array) {
	Array.getStatistics(array, min, max, mean, stdDev);
	return stdDev;
}


//Adds a scalar to all elements of an array
function addScalarToArray(array, scalar) {
	added_array=newArray(lengthOf(array));
	for (a=0; a<lengthOf(array); a++) {
		added_array[a]=array[a] + scalar;
	}
	return added_array;
}


//Get lower and upper percentile gray values
function get_percentile_values(lowerPercentile, upperPercentile, ignoreZeros) {
	resetMinAndMax();
	getRawStatistics(nPixels, mean, min, max, std, histogram);
	if(ignoreZeros == true) nPixels = nPixels - histogram[0];
	lowerTotal = 0;
	upperTotal = nPixels;

	i=0;
	if(ignoreZeros == true) i=1;
	while (lowerTotal < nPixels*lowerPercentile) {
		lowerTotal += histogram[i];
		//print("lower percentile: "+lowerTotal / (nPixels - histogram[0]));
		i++;
	}
	j=histogram.length-1;
	while (upperTotal > nPixels*upperPercentile) {
		upperTotal -= histogram[j];
		//print("upper: "+j+", "+upperTotal / (nPixels - histogram[0]));
		j--;
	}
	return newArray(i,j);
}
