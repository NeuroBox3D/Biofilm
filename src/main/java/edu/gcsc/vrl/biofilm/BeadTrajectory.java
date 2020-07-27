/* 
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS HEADER.
 *
 * Copyright (c) 2015–2020 G-CSC, Uni Frankfurt
 * 
 * This file is part of Visual Reflection Library (VRL).
 *
 * VRL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License version 3
 * as published by the Free Software Foundation.
 * 
 * see: http://opensource.org/licenses/LGPL-3.0
 *      file://path/to/VRL/src/eu/mihosoft/vrl/resources/license/lgplv3.txt
 *
 * VRL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * This version of VRL includes copyright notice and attribution requirements.
 * According to the LGPL this information must be displayed even if you modify
 * the source code of VRL. Neither the VRL Canvas attribution icon nor any
 * copyright statement/attribution may be removed.
 *
 * Attribution Requirements:
 *
 * If you create derived work you must do three things regarding copyright
 * notice and author attribution.
 *
 * First, the following text must be displayed on the Canvas:
 * "based on VRL source code". In this case the VRL canvas icon must be removed.
 * 
 * Second, the copyright notice must remain. It must be reproduced in any
 * program that uses VRL.
 *
 * Third, add an additional notice, stating that you modified VRL. A suitable
 * notice might read
 * "VRL source code modified by YourName 2012".
 * 
 * Note, that these requirements are in full accordance with the LGPL v3
 * (see 7. Additional Terms, b).
 *
 * Please cite the publication(s) listed below.
 *
 * Publications:
 *
 * M. Hoffer, C. Poliwoda, & G. Wittum. (2013). Visual reflection library:
 * a framework for declarative GUI programming on the Java platform.
 * Computing and Visualization in Science, 2013, 16(4),
 * 181–192. http://doi.org/10.1007/s00791-014-0230-y
 */
package edu.gcsc.vrl.biofilm;
import eu.mihosoft.vrl.annotation.ObjectInfo;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import javax.vecmath.Point4d;

/**
 * 
 * TODO testing TOL=1e-2
 *
 * @author kshitizmalhotra and Gillian Queisser
 */
@ObjectInfo(serializeParam=false)
public class BeadTrajectory {
    public int trajId;
    public File file;
    public double velocity;
    public List<Double> segmentLocalVelocities = new ArrayList<Double>();
    public double meanSegmentLocalVelocities;
    public double varianceSegmentLocalVelocities;
    public double density;
    public List<Double> segmentLocalDensities = new ArrayList<Double>();
    public double meanSegmentLocalDensities;
    public double varianceSegmentLocalDensities;
    public double weightedVelocity;
    public double localWeightedVelocity;
    public List<Double> segmentLocalWeightedVelocities = new ArrayList<Double>();
    public double meanSegmentLocalWeightedVelocities;
    public double varianceSegmentLocalWeightedVelocities;
    public int trajLength;
    public double physTrajLength;
    //public ImageVoxels voxels;  //? 
    public TrajBoundingBox box;
    public double boundingBoxVolume;
    public int point4dSize;
    public int maxSize;
    public int timestep;
    public List<Point4d> pointsOnTrajectory = new ArrayList<Point4d>(); // point in space-time (x,y,z,w=time)
    public List<Double> meanSquareDisplacement = new ArrayList<Double>();
    
    //public List<Point4d> trajects = new ArrayList<Point4d>();
            
   public BeadTrajectory() { 

    }
   
    public BeadTrajectory(double width){
        box.inputvoxelWidth = width;
    }
   
    public int getTrajId() {
        return trajId;
    }

    public void setTrajId(int trajId) {
        this.trajId = trajId;
    }
    public double getVelocity() {
        return velocity;
    }

    public void setVelocity(double velocity) {
        this.velocity = velocity;
    }
    
    public double getDensity() {
        return density;
    }
    public void setTimeStep(int timestep) {
        this.timestep = timestep;
    }

    public void setDensity(double density) {
        this.density = density;
    }
    
    public double getWeightedVelocity() {
        return weightedVelocity;
    }

    public void setWeightedVelocity(double weightedVelocity) {
        this.weightedVelocity = weightedVelocity;
    }
    
//    public double[] getSegmentWeightedLocalVelocity() {
//        return segmentLocalWeightedVelocities;
//    }

//    public void setSegmentWeightedLocalVelocity(double[] segmentLocalWeightedVelocities ) {
//        this.segmentLocalWeightedVelocities = segmentLocalWeightedVelocities;
//    }
    
    
     public double getLocalWeightedVelocity() {
        return localWeightedVelocity;
    }
    
    public void setLocalWeightedVelocity(double weightedVelocity) {
        
        this.localWeightedVelocity = localWeightedVelocity;
        
        throw new RuntimeException("Self asignment, does not work!");
    }
     
    public int getTrajLength() {
        return trajLength;
    }

    public void setTrajLength(int trajlength) {
        this.trajLength = trajlength;
    }
    public void setPoint4dSize(int size) {
        this.point4dSize = size;
    }
    
 //   public double getMeanSquareDisplacement(int trajId) {
 //       return meanSquareDisplacement;
 //   }

 //   public void setMeanSquareDisplacement(int trajId) {
 //       this.meanSquareDisplacement = meanSquareDisplacement;
 //   }
    
    public double getMeanSegmentLocalVelocities(int trajId){
        return meanSegmentLocalVelocities;
    }
    
    public void setMeanSegmentLocalVelocities(int trajId) {
        this.meanSegmentLocalVelocities = meanSegmentLocalVelocities;
        
        throw new RuntimeException("Self asignment, does not work!");
    }
    
    public double getVarianceSegmentLocalVelocities(int trajId){
        return varianceSegmentLocalVelocities;
    }
    
    public void setVarianceSegmentLocalVelocities(int trajId) {
        this.varianceSegmentLocalVelocities = varianceSegmentLocalVelocities;
        throw new RuntimeException("Self asignment, does not work!");
    }
    
    public double getMeanSegmentLocalDensities(int trajId){
        return meanSegmentLocalDensities;
    }
    
    public void setMeanSegmentLocalDensities(int trajId) {
        this.meanSegmentLocalDensities = meanSegmentLocalDensities;
        throw new RuntimeException("Self asignment, does not work!");
    }
    
    public double getVarianceSegmentLocalDensities(int trajId){
        return varianceSegmentLocalDensities;
    }
    
    public void setVarianceSegmentLocalDensities(int trajId) {
        this.varianceSegmentLocalDensities = varianceSegmentLocalDensities;
        throw new RuntimeException("Self asignment, does not work!");
    }
    
    public double getMeanSegmentLocalWeightedVelocities(int trajId){
        return meanSegmentLocalWeightedVelocities;
    }
    
    public void setMeanSegmentLocalWeightedVelocities(int trajId) {
        this.meanSegmentLocalWeightedVelocities = meanSegmentLocalWeightedVelocities;
        throw new RuntimeException("Self asignment, does not work!");
    }
    
    public double getVarianceSegmentLocalWeightedVelocities(int trajId){
        return varianceSegmentLocalWeightedVelocities;
    }
    
    public void setVarianceSegmentLocalWeightedVelocities(int trajId) {
        this.varianceSegmentLocalWeightedVelocities = varianceSegmentLocalWeightedVelocities;
        throw new RuntimeException("Self asignment, does not work!");
    }
    
    
     public  void computeVelocity(){
        double timeStep = this.timestep;
        double pixelWidth = box.getWidth();
        double pixelHeight = box.getHeight();
        double pixelDepth = box.getDepth();
        double velocity = 0.0;
        double totalDist = 0.0;
        int frameCurrent;
        double xCurrent;
        double yCurrent;
        double zCurrent;
        int frameNext;
        double xNext;
        double yNext;
        double zNext;
        
        
       // average velocity = length-of-segment / duration 
       for(int i = 0;i<this.pointsOnTrajectory.size()-1; i++){
            xCurrent = this.pointsOnTrajectory.get(i).x;
            yCurrent = this.pointsOnTrajectory.get(i).y;
            zCurrent = this.pointsOnTrajectory.get(i).z;
            frameNext = (int)this.pointsOnTrajectory.get(i+1).w;
            frameCurrent = (int)this.pointsOnTrajectory.get(i).w;
            xNext = this.pointsOnTrajectory.get(i+1).x;
            yNext = this.pointsOnTrajectory.get(i+1).y;
            zNext = this.pointsOnTrajectory.get(i+1).z;
            
              double dist = Math.sqrt((xCurrent-xNext)*(xCurrent-xNext)*pixelHeight*pixelHeight + (yCurrent-yNext)*(yCurrent-yNext)*pixelWidth*pixelWidth + (zCurrent-zNext)*(zCurrent-zNext)*pixelDepth*pixelDepth);
              double time = frameNext - frameCurrent;
              time = time * timeStep;
              velocity += dist/time;
              //fill list of velocities for each segment
              this.segmentLocalVelocities.add(dist/time);
              totalDist += dist;
            }
    
            // average = velocity/number of segments (= num points -1)
            velocity  = velocity/(this.getTrajLength()-1);

       this.velocity = velocity;
       this.physTrajLength = totalDist;
     }
     
    public void computeMeanSegmentLocalVelocity(){
        double sum = 0.0;
        //int l = Array.getLength(segmentLocalVelocities);
        for(int i = 0; i < this.pointsOnTrajectory.size()-1; i++){
            sum = sum + this.segmentLocalVelocities.get(i);
        }
        this.meanSegmentLocalVelocities = sum/(this.pointsOnTrajectory.size()-1);
    }
    
    public void computeVarianceSegmentLocalVelocity(){
        double sum = 0.0;
        //int l = Array.getLength(segmentLocalVelocities);
        for (int i = 0; i < this.pointsOnTrajectory.size()-1; i++){
            sum = sum + (this.segmentLocalVelocities.get(i) - this.meanSegmentLocalVelocities)*(this.segmentLocalVelocities.get(i) - this.meanSegmentLocalVelocities); 
        }
        this.varianceSegmentLocalVelocities = sum/(this.pointsOnTrajectory.size()-2);
    }
     
    public void computeTrajLength() {
            this.trajLength = pointsOnTrajectory.size();
        }
    public double getPhysTrajLength(){
        return physTrajLength; 
    }
    
    public void computeWeightedVelocity(){
        this.weightedVelocity = this.velocity * this.density;
    }

    public void computeLocalWeightedVelocity(ImageVoxels voxels){
       //computes the local density of a single segment of a trajectory and computes the weighted local velocity
        double timeStep = this.timestep;
        double pixelWidth = box.getWidth();
        double pixelHeight = box.getHeight();
        double pixelDepth = box.getDepth();
        double velocity = 0.0;              // TODO fix field hiding (velocity only used locally, no intention to set this.velocity)
        double localWeightedVelocity = 0.0; // TODO fix field hiding -> this.localWeightedVelocity is set at end of method
        //double totalDist = 0.0;
        int frameCurrent;
        double xCurrent;
        double yCurrent;
        double zCurrent;
        int frameNext;
        double xNext;
        double yNext;
        double zNext;
        double dist = 0.0;
        double time = 0.0;
        double localDensity = 0.0;
        
        // for each segment
        for(int i = 0;i<this.pointsOnTrajectory.size()-1; i++){
            xCurrent = this.pointsOnTrajectory.get(i).x;
            yCurrent = this.pointsOnTrajectory.get(i).y;
            zCurrent = this.pointsOnTrajectory.get(i).z;
            frameNext = (int)this.pointsOnTrajectory.get(i+1).w;
            frameCurrent = (int)this.pointsOnTrajectory.get(i).w;
            xNext = this.pointsOnTrajectory.get(i+1).x;
            yNext = this.pointsOnTrajectory.get(i+1).y;
            zNext = this.pointsOnTrajectory.get(i+1).z;
            
              dist = Math.sqrt((xCurrent-xNext)*(xCurrent-xNext)*pixelHeight*pixelHeight + (yCurrent-yNext)*(yCurrent-yNext)*pixelWidth*pixelWidth + (zCurrent-zNext)*(zCurrent-zNext)*pixelDepth*pixelDepth);
              time = frameNext - frameCurrent;
              time = time * timeStep;
              
              //compute the local density around the segment
              localDensity = computeLocalDensity(i, i+1, voxels);
              //fill list of local densities for each segment
              this.segmentLocalDensities.add(localDensity);
              
              localWeightedVelocity = dist/time * localDensity;
              
              //fill list of local weighted velocities for each segment
              this.segmentLocalWeightedVelocities.add(localWeightedVelocity);
              
              velocity += localWeightedVelocity;
              //totalDist += dist;
            }
    
            velocity  = velocity/(this.getTrajLength()-1);

       this.localWeightedVelocity = velocity;
      // this.physTrajLength = totalDist;
       
       
    }
    
    public void computeMeanSegmentLocalWeightedVelocity(){
        //double arrAvg = arr => arr.reduce((a,b) => a + b, 0)/arr.length;
        //double total = this.segmentLocalWeightedVelocities.reduce((previous,current) => current + previous);
        double sum = 0.0;
        //int l = Array.getLength(segmentLocalWeightedVelocities);
        for(int i =0; i < this.pointsOnTrajectory.size()-1; i++){
            sum = sum + this.segmentLocalWeightedVelocities.get(i);
        }
        this.meanSegmentLocalWeightedVelocities = sum/(this.pointsOnTrajectory.size()-1);
    }
    
    public void computeVarianceSegmentLocalWeightedVelocity(){
        double sum = 0.0;
        //int l = Array.getLength(segmentLocalWeightedVelocities);
        for(int i = 0; i < this.pointsOnTrajectory.size()-1; i++){
            sum = sum + (this.segmentLocalWeightedVelocities.get(i) - this.meanSegmentLocalWeightedVelocities)*(this.segmentLocalWeightedVelocities.get(i) - this.meanSegmentLocalWeightedVelocities);
        }
        this.varianceSegmentLocalWeightedVelocities = sum/(this.pointsOnTrajectory.size()-2);
    }
    
    public double computeLocalDensity(int point1, int point2, ImageVoxels voxels){
        //compute the local dimensions of the bounding box surrounding a single segment
        double x1 = pointsOnTrajectory.get(point1).x;
        double y1 = pointsOnTrajectory.get(point1).y;
        double z1 = pointsOnTrajectory.get(point1).z;
        
        double x2 = pointsOnTrajectory.get(point2).x;
        double y2 = pointsOnTrajectory.get(point2).y;
        double z2 = pointsOnTrajectory.get(point2).z;
        
        double xMin = 0.0;
        double yMin = 0.0;
        double zMin = 0.0;
        
        double xMax = 0.0;
        double yMax = 0.0;
        double zMax = 0.0;
        
        if (x1 < x2){
            xMin = x1;
            xMax = x2;
        } else{
            xMin = x2;
            xMax = x1;
        }
        
        if (y1 < y2){
            yMin = y1;
            yMax = y2;
        } else{
            yMin = y2;
            yMax = y1;
        }
        
        if (z1 < z2){
            zMin = z1;
            zMax = z2;
        } else{
            zMin = z2;
            zMax = z1;
        }
        
        double segmentWidth = java.lang.Math.abs(xMax - xMin);
        double segmentHeight = java.lang.Math.abs(yMax - yMin);
        double segmentDepth = java.lang.Math.abs(zMax - zMin);
        
        double segmentVolume = segmentWidth * segmentHeight * segmentDepth;
        
        //compute the local density
        double averageDensity = 0.0;
        
        int numVoxZ = (int)(segmentDepth);
                if(numVoxZ < 1){
                    numVoxZ = 1;
                }
                int numVoxY = (int)(segmentHeight);
                if(numVoxY < 1){
                    numVoxY = 1;
                }
                int numVoxX = (int)(segmentWidth);
                if(numVoxX < 1){
                    numVoxX = 1;
                }
                
        int numVoxel = (int)((numVoxX * numVoxY * numVoxZ)); 
        
         //  System.out.println("  numVox_z  "+numVox_z + "  numVox_y  "+numVox_y + "  numVox_x  "+numVox_x + "  z_min  "+z_min + "  y_min  "+y_min+"  x_min  "+x_min);
                    for (int zz = (int)zMin; zz < (int)zMin+numVoxZ; zz++) {
                        for (int yy = (int)yMin; yy < (int)yMin+numVoxY; yy++) {
                            for (int xx = (int)xMin; xx < (int)xMin+numVoxX; xx++) {
                                averageDensity  += voxels.getData()[zz][xx][yy];
                              // System.out.println("Density computed in computeDensity:    "+value);
                            }
                        }
                    }
                    
                    System.out.println("Number of voxels computed in computeLocalDensity:    "+numVoxel);
                    
                    //divide by number of voxels for average gray value
                    averageDensity /= numVoxel;
                     System.out.println("Density computed in computeLocalDensity:    " + averageDensity);
                    
                    //if necessary, rescale to [0,1] interval rather than [0,255] 
                    //value /= 255.0; // scale to [0,1], image has 8-bit (min=0,max=255)
                    
                    //on a 0-255 density scale, avoid density equal to zero, as weighted velocity will otherwise be infinity
                    if(averageDensity < 1){
                        averageDensity = 1.0;
                    }
        
        return averageDensity;
    }
    public void computeMeanSegmentLocalDensity(){
        double sum = 0.0;
        //int l = Array.getLength(segmentLocalDensities);
        
        // for each segment
        for(int i =0; i < this.pointsOnTrajectory.size()-1; i++){
            sum = sum + this.segmentLocalDensities.get(i);
        }
        this.meanSegmentLocalDensities = sum/(this.pointsOnTrajectory.size()-1);
    }
    
    /**
     * Variance of the segment level density. one traj for each particle. 
     */
    public void computeVarianceSegmentLocalDensity(){
        double sum = 0.0;
        //int l = Array.getLength(segmentLocalDensities);
        
        // for each segment
        for(int i = 0; i < this.pointsOnTrajectory.size()-1; i++){
            sum = sum + (this.segmentLocalDensities.get(i) - this.meanSegmentLocalDensities)*(this.segmentLocalDensities.get(i) - this.meanSegmentLocalDensities);
        }
        
        // why -2?
        // sigma^2 = sum (x_i - x_mean)^2 / (n-1), n = num segments
        this.varianceSegmentLocalDensities = sum/(this.pointsOnTrajectory.size()-2);
    }
    
    public void computeDimensions(){
        
        
        box.X_min = pointsOnTrajectory.get(0).x;
        box.X_max = pointsOnTrajectory.get(0).x;
        box.Y_min = pointsOnTrajectory.get(0).y;
        box.Y_max = pointsOnTrajectory.get(0).y;
        box.Z_min = pointsOnTrajectory.get(0).z;
        box.Z_max = pointsOnTrajectory.get(0).z;
            
            for(int i = 0;i<pointsOnTrajectory.size()-1;i++){
                
                if(this.pointsOnTrajectory.get(i+1).x>box.X_max){
                  box.X_max = this.pointsOnTrajectory.get(i+1).x;
                }else if(this.pointsOnTrajectory.get(i+1).x<box.X_min){
                  box.X_min = this.pointsOnTrajectory.get(i+1).x;
                }
                
                if(this.pointsOnTrajectory.get(i+1).y>box.Y_max){
                  box.Y_max = this.pointsOnTrajectory.get(i+1).y;
                }else if(this.pointsOnTrajectory.get(i+1).y<box.Y_min){
                  box.Y_min = this.pointsOnTrajectory.get(i+1).y;
                }
                    
                if(this.pointsOnTrajectory.get(i+1).z>box.Z_max){
                  box.Z_max = this.pointsOnTrajectory.get(i+1).z;
                }else if(this.pointsOnTrajectory.get(i+1).z<box.Z_min){
                  box.Z_min = this.pointsOnTrajectory.get(i+1).z;
                }
            }
            
            box.cubeWidth = box.X_max - box.X_min;
            box.cubeHeight = box.Y_max - box.Y_min;
            box.cubeDepth = box.Z_max - box.Z_min;
            this.boundingBoxVolume = box.cubeWidth * box.cubeHeight * box.cubeDepth;
        }
    public void computeDensity(ImageVoxels voxels){
           
        //compute the number of voxels in bounding box
        double value = 0;
                int numVox_z = (int)(box.cubeDepth);
                if(numVox_z < 1){
                    numVox_z = 1;
                }
                int numVox_y = (int)(box.cubeHeight);
                if(numVox_y < 1){
                    numVox_y = 1;
                }
                int numVox_x = (int)(box.cubeWidth);
                if(numVox_x < 1){
                    numVox_x = 1;
                }
                int numVoxel = (int)((numVox_x * numVox_y * numVox_z)) ; 
               
        //compute the density in bounding box. Iterate from x_min, y_min, z_min to x_max, y_max, z_max and sum up the gray values        
                int z_min = (int)(box.Z_min); 
                int y_min = (int)(box.Y_min); 
                int x_min = (int)(box.X_min); 
                
                  System.out.println("  numVox_z  "+numVox_z + "  numVox_y  "+numVox_y + "  numVox_x  "+numVox_x + "  z_min  "+z_min + "  y_min  "+y_min+"  x_min  "+x_min);
                    for (int zz = z_min; zz < z_min+numVox_z; zz++) {
                        for (int yy = y_min; yy < y_min+numVox_y; yy++) {
                            for (int xx = x_min; xx < x_min+numVox_x; xx++) {
                                value  += voxels.getData()[zz][xx][yy];
                               System.out.println("Density computed in computeDensity:    "+value);
                            }
                        }
                    }
                    
                    System.out.println("Number of voxels computed in computeDensity:    "+numVoxel);
                    
                    //divide by number of voxels for average gray value
                    value /= numVoxel;
                     System.out.println("Density computed in computeDensity:    "+value);
                    
                    //if necessary, rescale to [0,1] interval rather than [0,255] 
                    //value /= 255.0; // scale to [0,1], image has 8-bit (min=0,max=255)
                    
                    //on a 0-255 density scale, avoid density equal to zero, as weighted velocity will otherwise be infinity
                    if(value < 1){
                        value = 1;
                    }
                    
                    this.density = value;
    }
    
    /**
     * Used to determine brownian motion vs. active transport. For each point.
     */
    public void computeMeanSquareDisplacement(){
        //compute the mean square displacement MSD = 1/N * \sum_i{(x0 - xi)^2}
        
        double x0 = this.pointsOnTrajectory.get(0).x;
        double y0 = this.pointsOnTrajectory.get(0).y;
        double z0 = this.pointsOnTrajectory.get(0).z;
        
        double xCurrent = 0.0;
        double yCurrent = 0.0;
        double zCurrent = 0.0;
        
        double distance = 0.0;
        double MSD = 0.0;
        
        // for each trajectory point
        for(int i = 0; i<pointsOnTrajectory.size(); i++){
            xCurrent = this.pointsOnTrajectory.get(i).x;
            yCurrent = this.pointsOnTrajectory.get(i).y;
            zCurrent = this.pointsOnTrajectory.get(i).z;
            
            distance = Math.sqrt((xCurrent-x0)*(xCurrent-x0) + (yCurrent-y0)*(yCurrent-y0) + (zCurrent-z0)*(zCurrent-z0));
            MSD += distance * distance;
            this.meanSquareDisplacement.add(MSD/(i+1));
        }
        
        //meanSquareDisplacement /= this.trajLength;
    }
}
