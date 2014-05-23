//
//  main.m
//  MomentTensorModelStatistics
//
//  Created by Jeffrey J. Early on 5/22/14.
//  Copyright (c) 2014 Jeffrey J. Early. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>
#import <GLNumericalModelingKit/GLOperationOptimizer.h>
#import "MomentTensorModels.h"

void displayKappaSimple( GLFunction *t, GLFunction *xPosition, GLFunction *yPosition, GLFloat kappa)
{
	GLFunction *meanSquareSeparation = [[[xPosition times: xPosition] plus: [yPosition times: yPosition]] mean: 1];
	
	GLFloat a = meanSquareSeparation.pointerValue[0];
	GLFloat b = meanSquareSeparation.pointerValue[meanSquareSeparation.nDataPoints-1];
	
	GLFloat kappaDeduced = (0.25)*(b-a)/t.pointerValue[t.nDataPoints-1];
	NSLog(@"kappa: %f, actual kappa: %f", kappa, kappaDeduced);
}

int main(int argc, const char * argv[])
{

	@autoreleasepool {
		GLFloat floatSpacing = 10;
		GLFloat maxTime = 6*86400;
		GLFloat timeStep = 30*60;
		NSInteger nParticles = 10;
	    GLEquation *equation = [[GLEquation alloc] init];
		GLDimension *floatDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: nParticles domainMin: 1 length: nParticles];
		GLFunction *xPosition = [GLFunction functionOfRealTypeWithDimensions: @[floatDim] forEquation: equation];
		GLFunction *yPosition = [GLFunction functionOfRealTypeWithDimensions: @[floatDim] forEquation: equation];
		
		// Layout the drifters in a cross pattern, just the real drifters
		NSUInteger iFloat = 0;
		GLFloat length = ((GLFloat) nParticles/2 - 1)*floatSpacing;
		for (NSInteger i=0; i<nParticles/2; i++) {
			xPosition.pointerValue[iFloat] = ((GLFloat) i)*floatSpacing - length/2;
			yPosition.pointerValue[iFloat] = ((GLFloat) i)*0;
			iFloat++;
		}
		for (NSInteger i=0; i<nParticles/2; i++) {
			//if (i==0) continue;
			xPosition.pointerValue[iFloat] = ((GLFloat) i)*0;
			yPosition.pointerValue[iFloat] = ((GLFloat) i)*floatSpacing - length/2;
			iFloat++;
		}
	    
		GLFloat kappa = 0.2; // m^2/s
        GLFloat norm = sqrt(timeStep*2*kappa);
        norm = sqrt(36./10.)*norm/timeStep; // the integrator multiplies by deltaT, so we account for that here.
        // RK4: dt/3 f(0) + dt/6 f(1) + dt/6 *f(4) + dt/3*f(3)
        // sqrt of ( (1/3)^2 + (1/6)^ + (1/6)^2 + (1/3)^2 )
        
        NSArray *y=@[xPosition, yPosition];
        GLRungeKuttaOperation *integrator = [GLRungeKuttaOperation rungeKutta4AdvanceY: y stepSize: timeStep fFromTY:^(GLScalar *time, NSArray *yNew) {
            GLFunction *xStep = [GLFunction functionWithNormallyDistributedValueWithDimensions: @[floatDim] forEquation: equation];
            GLFunction *yStep = [GLFunction functionWithNormallyDistributedValueWithDimensions: @[floatDim] forEquation: equation];
            xStep = [xStep times: @(norm)];
            yStep = [yStep times: @(norm)];
            
			return @[xStep, yStep];
		}];
		
		GLDimension *tDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: 1+round(maxTime/timeStep) domainMin: 0 length: maxTime];
		GLFunction *t = [GLFunction functionOfRealTypeFromDimension: tDim withDimensions: @[tDim] forEquation: equation];
		
		for (NSUInteger i=0; i<10; i++) {
			
			NSArray *newPositions = [integrator integrateAlongDimension: tDim];
					
			displayKappaSimple( t, newPositions[0],  newPositions[1], kappa);
			
			MomentTensorModels *models = [[MomentTensorModels alloc] initWithXPositions: newPositions[0] yPositions:newPositions[1] time: t];
			NSArray *result = [models bestFitToDiffusivityModel];
			
			GLScalar *minError = result[0];
			GLScalar *minKappa = result[1];
			
			NSLog(@"diffusivity model\t\t\t\t\terror: %f (kappa)=(%.4f)", *(minError.pointerValue), *(minKappa.pointerValue));
			
			result = [models bestFitToStrainDiffusivityModel];
			
			minError = result[0];
			minKappa = result[1];
			GLScalar *minSigma = result[2];
			GLScalar *minTheta = result[3];
			
			NSLog(@"strain-diffusivity model\t\t\t\terror: %f (kappa,sigma,theta)=(%.4f,%.3g,%.1f)", *(minError.pointerValue), *(minKappa.pointerValue),*(minSigma.pointerValue),(*(minTheta.pointerValue))*180./M_PI);
            
            
			result = [models bestFitToVorticityStrainDiffusivityModel];
			
			minError = result[0];
			minKappa = result[1];
			minSigma = result[2];
			minTheta = result[3];
			GLScalar *minZeta = result[4];
						
			NSLog(@"strain-vorticity-diffusivity model\terror: %f (kappa,sigma,theta)=(%.4f,%.3g,%.1f,%.3g)", *(minError.pointerValue), *(minKappa.pointerValue),*(minSigma.pointerValue),(*(minTheta.pointerValue))*180./M_PI, *(minZeta.pointerValue));

		}
	}
    return 0;
}



