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


/**
 *
 * @author kshitizmalhotra
 */
public class TrajBoundingBox extends BeadTrajectory{
    public double cubeWidth;
    public ImageVoxels voxels;
    public double cubeHeight;
    public double cubeDepth;
    public double inputvoxelWidth;
    public double inputvoxelHeight;
    public double inputvoxelDepth;
    public double X_min;
    public double Y_min;
    public double Z_min;
    public double X_max;
    public double Y_max;
    public double Z_max;
    
    

public TrajBoundingBox(double cubeWidth, double cubeHeight, double  cubeDepth) {
        this.cubeWidth = cubeWidth;
        this.cubeHeight = cubeHeight;
        this.cubeDepth = cubeDepth;
        
    }


public TrajBoundingBox(){
    //this.cubeWidth = 0.0;
    //this.cubeHeight = 0.0;
    //this.cubeDepth = 0.0;
    //this.x_min = 0.0;
    //this.y_min = 0.0;
    //this.z_min = 0.0;
    //this.inputvoxelWidth = 0.0;
    //this.inputvoxelHeight = 0.0;
    //this.inputvoxelDepth = 0.0;
}
    public double getcubeWidth() {
        return cubeWidth;
    }

    public void setcubeWidth(double cubeWidth) {
        this.cubeWidth = cubeWidth;
    }
    public double getcubeHeight() {
        return cubeHeight;
    }

    public void setcubeHeight(double  cubeHeight) {
        this.cubeHeight  = cubeHeight;
    }   
    public double  getcubeDepth() {
        return cubeDepth;
    }

    public void setcubeDepth(double cubeDepth) {
        this.cubeDepth = cubeDepth;
    }
    public void setx_min(double x_min) {
        this.X_min = x_min;
    }
    public void sety_min(double y_min) {
        this.Y_min = y_min;
    }
    public void setz_min(double z_min) {
        this.Z_min = z_min;
    }
    public double getWidth() {
        return inputvoxelWidth;
    }

    public void setWidth(double width_ui) {
        this.inputvoxelWidth = width_ui;
    }

    public double getHeight() {
        return inputvoxelHeight;
    }

    public void setHeight(double height_ui) {
        this.inputvoxelHeight = height_ui;
    }

    public double  getDepth() {
        return inputvoxelDepth;
    }

    public void setDepth(double depth_ui) {
        this.inputvoxelDepth = depth_ui;
    }
    /*public void computeDimensions(){
        
        
        X_min = b.pointsOnTrajectory.get(0).x;
        System.out.print("X_min is "+X_min);
        X_max = b.pointsOnTrajectory.get(0).x;
        Y_min = b.pointsOnTrajectory.get(0).y;
        Y_max = b.pointsOnTrajectory.get(0).y;

        Z_min = b.pointsOnTrajectory.get(0).z;
        Z_max = b.pointsOnTrajectory.get(0).z;
            
            for(int i = 0;i<b.pointsOnTrajectory.size()-1;i++){
                
                if(b.pointsOnTrajectory.get(i+1).x>X_max){
                  X_max = b.pointsOnTrajectory.get(i+1).x;
                }else if(b.pointsOnTrajectory.get(i+1).x<X_min){
                  X_min = b.pointsOnTrajectory.get(i+1).x;
                }
                
                if(this.b.pointsOnTrajectory.get(i+1).y>Y_max){
                  Y_max = b.pointsOnTrajectory.get(i+1).y;
                }else if(b.pointsOnTrajectory.get(i+1).y<Y_min){
                  Y_min = b.pointsOnTrajectory.get(i+1).y;
                }
                    
                if(this.b.pointsOnTrajectory.get(i+1).z>Z_max){
                  Z_max = b.pointsOnTrajectory.get(i+1).z;
                }else if(b.pointsOnTrajectory.get(i+1).z<Z_min){
                  Z_min = b.pointsOnTrajectory.get(i+1).z;
                }
            }
            
            this.cubeWidth = X_max - X_min;
            this.cubeHeight = Y_max - Y_min;
            this.cubeDepth = Z_max - Z_min;
            
        }*/
    
    
   /* public void computeDensity(){
           double value = 0;
           int numVoxel = (int)((this.cubeDepth*this.cubeHeight*this.cubeWidth)/(this.getWidth()*this.getDepth()*this.getHeight())) ;
                int numVox_z = (int)(this.cubeDepth/this.inputvoxelDepth);
                int numVox_y = (int)(this.cubeHeight/this.inputvoxelHeight);
                int numVox_x = (int)(this.cubeWidth/this.inputvoxelWidth);
                
                    for (int zz = 0; zz < numVox_z; zz++) {
                        for (int yy = 0; yy < numVox_y; yy++) {
                            for (int xx = 0; xx < numVox_x; xx++) {
                                value  += b.voxels.getData()[zz][xx][yy];
                            }
                        }
                    }
                    value /= numVoxel;
                    value /= 255.0; // scale to [0,1], image has 8-bit (min=0,max=255)
                    b.density = value;
    }

}*/



///    OR THIS 

 /*public void computeDensity(){
           double value = 0;
           int numVoxel = (int)((this.cubeDepth*this.cubeHeight*this.cubeWidth)/(this.getWidth()*this.getDepth()*this.getHeight())) ;
                double numVox_z = this.cubeDepth/this.inputvoxelDepth;
                double numVox_y = this.cubeHeight/this.inputvoxelHeight;
                double numVox_x = this.cubeWidth/this.inputvoxelWidth;
                
                
                  /* System.out.println("  something here  ");
                    for (double zz = Z_min; zz < numVox_z; zz = zz + this.inputvoxelDepth) {
                        System.out.println("  something here  ");
                        for (double yy = Y_min; yy < numVox_y; yy = yy + this.inputvoxelHeight) {
                            System.out.println("  something here  ");
                            for (double xx = X_min; xx < numVox_x; xx = xx + this.inputvoxelWidth) {
                                System.out.println("  something here  ");
                                System.out.println("  something else  ");
                                value  += voxels.getData()[(int)zz][(int)xx][(int)yy];
                                System.out.println("    "+voxels.getData()[(int)zz][(int)xx][(int)yy]);
                            }
                        }
                    }
                  System.out.println("  something here  ");
      
                                System.out.println("    "+voxels.getData()[0][0][0]);
                            
                  
                    value /= numVoxel;
                    value /= 255.0; // scale to [0,1], image has 8-bit (min=0,max=255)
                    density = value;
    }*/

}