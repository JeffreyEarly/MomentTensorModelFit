//
//  main.m
//  MomentTensorFit
//
//  Created by Jeffrey Early on 9/16/13.
//
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>
#import <GLNumericalModelingKit/GLOperationOptimizer.h>
#import "MomentTensorModels.h"
#import "MomentEllipseFileReader.h"
#import "EllipseErrorOperation.h"

int main(int argc, const char * argv[])
{
	
    @autoreleasepool {
        
        GLEquation *equation = [[GLEquation alloc] init];
        
		
		NSFileManager *fileManager = [[NSFileManager alloc] init];
		NSString *folderPath = @"/Users/jearly/Documents/LatMix/drifters/ObservationalData/griddedRhoDrifterMomementEllipses/";
        folderPath = @"/Users/jearly/Documents/LatMix/drifters/synthetic/moment-ellipses/synthetic-diffusive/";
        folderPath = @"/Users/jearly/Documents/LatMix/drifters/synthetic/moment-ellipses/synthetic-strained-diffusive/";
		folderPath = @"/Users/jearly/Documents/LatMix/drifters/synthetic/moment-ellipses/synthetic-vorticity-strained-diffusive/";
		NSArray *ellipseFiles = [fileManager contentsOfDirectoryAtPath: folderPath error: nil];
		
		NSMutableString *outputData = [NSMutableString stringWithFormat: @""];
		NSUInteger i=1;
		
		for ( NSString *filename in ellipseFiles)
		{
			if ( ![filename.pathExtension isEqualToString: @"txt"]) continue;
			MomentEllipseFileReader *file = [[MomentEllipseFileReader alloc] initWithURL: [NSURL fileURLWithPath: [NSString stringWithFormat: @"%@%@", folderPath, filename]] equation: equation];
			if (!file) continue;
			
			NSArray *drifterIDs = [[filename stringByDeletingPathExtension] componentsSeparatedByString: @"_"];
			NSMutableString *idString = [NSMutableString stringWithFormat: @"drifterIDs{%lu}=[", i];
			for (NSString *name in drifterIDs) {
				if ([drifterIDs indexOfObject: name] == 0) {
					continue; // the first string is just 'drifter'.
				}
				else if ([drifterIDs indexOfObject: name] == drifterIDs.count-1) {
					[idString appendFormat: @"%@", name];
				} else {
					[idString appendFormat: @"%@;", name];
				}
			}
			[idString appendString: @"]; "];
			[outputData appendString: idString];
			
			
			GLFloat kappaScale = 0.1;
			GLScalar *kappa = [GLScalar scalarWithValue: log(0.1/kappaScale) forEquation: equation];
			GLScalar *kappaDelta = [GLScalar scalarWithValue: 2.0 forEquation: equation];
			
			GLMinimizationOperation *minimizer = [[GLMinimizationOperation alloc] initAtPoint: @[kappa] withDeltas: @[kappaDelta] forFunction:^(NSArray *xArray) {
				GLScalar *kappaUnscaled = [[xArray[0] exponentiate] times: @(kappaScale)];
				NSArray *tensorComps = diffusivityModel( file.Mxx0, file.Myy0, file.Mxy0, file.t, kappaUnscaled);
				NSArray *ellipseComps = tensorCompsToEllipseComps( tensorComps );
		 
				EllipseErrorOperation *error = [[EllipseErrorOperation alloc] initWithParametersFromEllipseA: ellipseComps ellipseB:@[file.a, file.b, file.angle]];
			
				return error.result[0];
			}];
			
			GLScalar *minKappa = [[minimizer.result[0] exponentiate] times: @(kappaScale)];
			GLScalar *minError = minimizer.result[1];
			
			[outputData appendFormat: @"model1_error(%lu)=%g; model1_kappa(%lu)=%g; ", i, *(minError.pointerValue), i, *(minKappa.pointerValue)];
			
			NSLog(@"%@---diffusivity model total error: %f @ (kappa)=(%.4f)", filename.lastPathComponent, *(minError.pointerValue), *(minKappa.pointerValue));
			
            
            
            
			// kappaDelta is carefully chosen. It needs to represent a sort of 'size of parameter space' that we want to explore.
			// So here, we make it move around in fairly big chunks, like orders of magnitude.
			// Yup, big steps are the best, for ALL parameters.
			kappa = [GLScalar scalarWithValue: log(.1/kappaScale) forEquation: equation];
			kappaDelta = [GLScalar scalarWithValue: 3.0 forEquation: equation];
			
			GLFloat sigmaScale = 4.0E-6;
			GLScalar *sigma = [GLScalar scalarWithValue: log(1E-6/sigmaScale) forEquation: equation];
			GLScalar *sigmaDelta = [GLScalar scalarWithValue: 0.5 forEquation: equation];
			
			GLScalar *theta = [GLScalar scalarWithValue: 0.0*M_PI/180. forEquation: equation];
			GLScalar *thetaDelta = [GLScalar scalarWithValue: 45.*M_PI/180. forEquation: equation];
			
			minimizer = [[GLMinimizationOperation alloc] initAtPoint: @[kappa, sigma, theta] withDeltas: @[kappaDelta, sigmaDelta, thetaDelta] forFunction:^(NSArray *xArray) {
				GLScalar *kappaUnscaled = [[xArray[0] exponentiate] times: @(kappaScale)];
				GLScalar *sigmaUnscaled = [[xArray[1] exponentiate] times: @(sigmaScale)];
				NSArray *tensorComps = strainDiffusivityModel( file.Mxx0, file.Myy0, file.Mxy0, file.t, kappaUnscaled, sigmaUnscaled, xArray[2]);
				NSArray *ellipseComps = tensorCompsToEllipseComps( tensorComps );
				
				EllipseErrorOperation *error = [[EllipseErrorOperation alloc] initWithParametersFromEllipseA: ellipseComps ellipseB:@[file.a, file.b, file.angle]];
				
				return error.result[0];
			}];
			
            NSArray *results = minimizer.result;
            
			minKappa = [[minimizer.result[0] exponentiate] times: @(kappaScale)];
			GLScalar *minSigma = [[minimizer.result[1] exponentiate] times: @(sigmaScale)];
			GLScalar *minTheta = minimizer.result[2];
			minError = minimizer.result[3];
			
			[outputData appendFormat: @"model2_error(%lu)=%g; model2_kappa(%lu)=%g; model2_sigma(%lu)=%g; model2_theta(%lu)=%g;\n", i, *(minError.pointerValue), i, *(minKappa.pointerValue), i, *(minSigma.pointerValue), i, *(minTheta.pointerValue)];
			
			NSLog(@"%@---strain-diffusivity model total error: %f (kappa,sigma,theta)=(%.4f,%.3g,%.1f)", filename.lastPathComponent, *(minError.pointerValue), *(minKappa.pointerValue),*(minSigma.pointerValue),(*(minTheta.pointerValue))*180./M_PI);
            
            
            sigmaScale = 1e-5;
            sigma = [GLScalar scalarWithValue: log(1E-5/sigmaScale) forEquation: equation];
            
			
			
            GLFloat sScale = 1E-6;
			GLScalar *s = [GLScalar scalarWithValue: log(1e-6/sScale) forEquation: equation];
			GLScalar *sDelta = [GLScalar scalarWithValue: 0.5 forEquation: equation];
			
			GLScalar *alpha = [GLScalar scalarWithValue: 0 forEquation: equation];
			GLScalar *alphaDelta = [GLScalar scalarWithValue: 1.0 forEquation: equation];
            
            // initializing with the results from the previous calculation. we can use log searches if we look for both positive and negative vorticity.
            GLMinimizationOperation *minimizer_positive = [[GLMinimizationOperation alloc] initAtPoint: @[kappa, s, theta, alpha] withDeltas: @[kappaDelta, sDelta, thetaDelta, alphaDelta] forFunction:^(NSArray *xArray) {
				GLScalar *kappaUnscaled = [[xArray[0] exponentiate] times: @(kappaScale)];
				GLScalar *sUnscaled = [[xArray[1] exponentiate] times: @(sScale)];
				GLScalar *sigmaUnscaled = [sUnscaled times: [xArray[3] cosh]];
                GLScalar *zetaUnscaled = [sUnscaled times: [xArray[3] sinh]];
				NSArray *tensorComps = strainVorticityDiffusivityModel( file.Mxx0, file.Myy0, file.Mxy0, file.t, kappaUnscaled, sigmaUnscaled, xArray[2], zetaUnscaled);
				NSArray *ellipseComps = tensorCompsToEllipseComps( tensorComps );
				
				EllipseErrorOperation *error = [[EllipseErrorOperation alloc] initWithParametersFromEllipseA: ellipseComps ellipseB:@[file.a, file.b, file.angle]];
				
				return error.result[0];
			}];
            NSArray *posResults = minimizer_positive.result;
            
			alphaDelta = [GLScalar scalarWithValue: -1.0 forEquation: equation];
			GLMinimizationOperation *minimizer_negative = [[GLMinimizationOperation alloc] initAtPoint: @[kappa, s, theta, alpha] withDeltas: @[kappaDelta, sDelta, thetaDelta, alphaDelta] forFunction:^(NSArray *xArray) {
				GLScalar *kappaUnscaled = [[xArray[0] exponentiate] times: @(kappaScale)];
				GLScalar *sUnscaled = [[xArray[1] exponentiate] times: @(sScale)];
				GLScalar *sigmaUnscaled = [sUnscaled times: [xArray[3] cosh]];
                GLScalar *zetaUnscaled = [sUnscaled times: [xArray[3] sinh]];
				NSArray *tensorComps = strainVorticityDiffusivityModel( file.Mxx0, file.Myy0, file.Mxy0, file.t, kappaUnscaled, sigmaUnscaled, xArray[2], zetaUnscaled);
				NSArray *ellipseComps = tensorCompsToEllipseComps( tensorComps );
				
				EllipseErrorOperation *error = [[EllipseErrorOperation alloc] initWithParametersFromEllipseA: ellipseComps ellipseB:@[file.a, file.b, file.angle]];
				
				return error.result[0];
			}];
            NSArray *negResults = minimizer_negative.result;
            
			GLScalar *minPos = posResults[4];
			GLScalar *minNeg = negResults[4];
			GLFloat pos = *(minPos.pointerValue);
			GLFloat neg = *(minNeg.pointerValue);
			
			results = pos < neg ? posResults : negResults;
			
			minKappa = [[results[0] exponentiate] times: @(kappaScale)];
			GLScalar *minS = [[results[1] exponentiate] times: @(sScale)];
			minSigma = [minS times: [results[3] cosh]];
            minTheta = results[2];
            GLScalar *minZeta = [minS times: [results[3] sinh]];
			minError = results[4];
			
			[outputData appendFormat: @"model3_error(%lu)=%g; model3_kappa(%lu)=%g; model3_sigma(%lu)=%g; model3_theta(%lu)=%g; model3_zeta(%lu)=%g;\n", i, *(minError.pointerValue), i, *(minKappa.pointerValue), i, *(minSigma.pointerValue), i, *(minTheta.pointerValue), i, *(minZeta.pointerValue)];
			
			NSLog(@"%@---strain-vorticity-diffusivity model total error: %f (kappa,sigma,theta)=(%.4f,%.3g,%.1f,%.3g)", filename.lastPathComponent, *(minError.pointerValue), *(minKappa.pointerValue),*(minSigma.pointerValue),(*(minTheta.pointerValue))*180./M_PI, *(minZeta.pointerValue));
            
			i++;
        }
		
		[outputData writeToFile: [NSString stringWithFormat: @"%@BestFitParameters_area_divergence.m", folderPath] atomically: YES encoding: NSUTF8StringEncoding error: nil];
    }
    return 0;
}
