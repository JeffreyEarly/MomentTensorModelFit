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
		GLFloat floatSpacing = 500;
		GLFloat maxTime = 6*86400;
		GLFloat timeStep = 30*60;
		NSInteger nParticles = 10;
	    GLEquation *equation = [[GLEquation alloc] init];
		GLDimension *floatDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: nParticles domainMin: 1 length: nParticles];
		floatDim.name = @"floats";
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
		
		NSMutableString *outputData = [NSMutableString string];
		for (NSUInteger iModel=0; iModel<3; iModel++)
		{
			GLFloat kappa = 0.2; // m^2/s
			
			NSString *name;
			NSUInteger totalIterations = 500;
			NSArray * (^addUV) (GLFunction *,GLFunction *, GLFunction *,GLFunction *);
			if (iModel == 0) {
				kappa = 0.569380;
				name = @"syntheticDiffusive";
				addUV = ^( GLFunction *xpos, GLFunction *ypos, GLFunction *u,GLFunction *v ) {
					return @[u,v];
				};
				[outputData appendFormat:@"%@.truth = struct('kappa', %g, 'sigma', 0.0, 'theta', 0.0, 'zeta', 0.0, 'error', 0.0);\n", name, kappa];
			} else if (iModel == 1) {
				name = @"syntheticStrainDiffusive";
				kappa = 0.205619;
				GLFloat sigma = 3.49415e-6;
				GLFloat theta = -32.359313*M_PI/180.;
				GLFloat sigma_n = sigma*cos(2.*theta);
				GLFloat sigma_s = sigma*sin(2.*theta);
				addUV = ^( GLFunction *xpos, GLFunction *ypos, GLFunction *u,GLFunction *v ) {
					GLFunction *u2 = [[xpos times: @(sigma_n/2.)] plus: [ypos times: @(sigma_s/2.)]];
					GLFunction *v2 = [[xpos times: @(sigma_s/2.)] plus: [ypos times: @(-sigma_n/2.)]];
					return @[[u plus: u2],[v plus: v2]];
				};
				[outputData appendFormat:@"%@.truth = struct('kappa', %g, 'sigma', %g, 'theta', %g, 'zeta', 0.0, 'error', 0.0);\n", name, kappa, sigma, theta];
			} else {
				name = @"syntheticVorticityStrainDiffusive";
//				GLFloat s = 8e-6;
//				GLFloat alpha = 0.75;
//				GLFloat sigma = s*cosh(alpha);
//				GLFloat zeta = s*sinh(alpha);
				kappa = 0.201481;
				GLFloat sigma = 3.50241e-06;
				GLFloat zeta = -3.75091e-08;
				GLFloat theta = -32.278997*M_PI/180.;
				GLFloat sigma_n = sigma*cos(2.*theta);
				GLFloat sigma_s = sigma*sin(2.*theta);
				
				addUV = ^( GLFunction *xpos, GLFunction *ypos, GLFunction *u,GLFunction *v ) {
					GLFunction *u2 = [[xpos times: @(sigma_n/2.)] plus: [ypos times: @((sigma_s-zeta)/2.)]];
					GLFunction *v2 = [[xpos times: @((sigma_s + zeta)/2.)] plus: [ypos times: @(-sigma_n/2.)]];
					return @[[u plus: u2],[v plus: v2]];
				};
				[outputData appendFormat:@"%@.truth = struct('kappa', %g, 'sigma', %g, 'theta', %g, 'zeta', %g, 'error', 0.0);\n", name, kappa, sigma, theta, zeta];
			}
			
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
				
				return addUV(y[0],y[1],xStep,yStep);
			}];
			
			GLDimension *tDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: 1+round(maxTime/timeStep) domainMin: 0 length: maxTime];
			tDim.name = @"t";
			GLFunction *t = [GLFunction functionOfRealTypeFromDimension: tDim withDimensions: @[tDim] forEquation: equation];
			
			
			[outputData appendFormat:@"%@.model1(%lu) = struct('kappa', 0.0, 'sigma', 0.0, 'theta', 0.0, 'zeta', 0.0, 'error', 0.0);\n", name, totalIterations];
			[outputData appendFormat:@"%@.model2(%lu) = struct('kappa', 0.0, 'sigma', 0.0, 'theta', 0.0, 'zeta', 0.0, 'error', 0.0);\n", name, totalIterations];
			[outputData appendFormat:@"%@.model3(%lu) = struct('kappa', 0.0, 'sigma', 0.0, 'theta', 0.0, 'zeta', 0.0, 'error', 0.0);\n", name, totalIterations];
			for (NSUInteger i=0; i<totalIterations; i++) {
				
				NSArray *newPositions = [integrator integrateAlongDimension: tDim];
						
				displayKappaSimple( t, newPositions[0],  newPositions[1], kappa);
				
				MomentTensorModels *models = [[MomentTensorModels alloc] initWithXPositions: newPositions[0] yPositions:newPositions[1] time: t];
				NSArray *result = [models bestFitToDiffusivityModel];
				
				GLScalar *minError = result[0];
				GLScalar *minKappa = result[1];
				
				[outputData appendFormat: @"%@.model1(%lu) = struct('kappa', %g, 'sigma', 0.0, 'theta', 0.0, 'zeta', 0.0, 'error', %g);\n", name, i+1, *(minKappa.pointerValue),*(minError.pointerValue)];
				NSLog(@"diffusivity model\t\t\t\t\terror: %f (kappa)=(%.4f)", *(minError.pointerValue), *(minKappa.pointerValue));
				
				result = [models bestFitToStrainDiffusivityModel];
				
				minError = result[0];
				minKappa = result[1];
				GLScalar *minSigma = result[2];
				GLScalar *minTheta = result[3];
				
				[outputData appendFormat: @"%@.model2(%lu) = struct('kappa', %g, 'sigma', %g, 'theta', %g, 'zeta', 0.0, 'error', %g);\n", name, i+1, *(minKappa.pointerValue),*(minSigma.pointerValue),*(minTheta.pointerValue),*(minError.pointerValue)];
				NSLog(@"strain-diffusivity model\t\t\t\terror: %f (kappa,sigma,theta)=(%.4f,%.3g,%.1f)", *(minError.pointerValue), *(minKappa.pointerValue),*(minSigma.pointerValue),(*(minTheta.pointerValue))*180./M_PI);
				
				
				result = [models bestFitToVorticityStrainDiffusivityModelWithStartPoint: @[minKappa, minSigma, minTheta]];
				
				minError = result[0];
				minKappa = result[1];
				minSigma = result[2];
				minTheta = result[3];
				GLScalar *minZeta = result[4];
							
				[outputData appendFormat: @"%@.model3(%lu) = struct('kappa', %g, 'sigma', %g, 'theta', %g, 'zeta', %g, 'error', %g);\n", name, i+1, *(minKappa.pointerValue),*(minSigma.pointerValue),*(minTheta.pointerValue),*(minZeta.pointerValue),*(minError.pointerValue)];
				NSLog(@"strain-vorticity-diffusivity model\terror: %f (kappa,sigma,theta)=(%.4f,%.3g,%.1f,%.3g)", *(minError.pointerValue), *(minKappa.pointerValue),*(minSigma.pointerValue),(*(minTheta.pointerValue))*180./M_PI, *(minZeta.pointerValue));

			}
		}
		[outputData writeToFile: @"/Users/jearly/Documents/LatMix/drifters/synthetic/BestFitEllipseAreaDivergenceTotalSum.m" atomically: YES encoding: NSUTF8StringEncoding error: nil];
	}
    return 0;
}



