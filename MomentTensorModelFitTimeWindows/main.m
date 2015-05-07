//
//  main.m
//  MomentTensorModelFitTimeWindows
//
//  Created by Jeffrey J. Early on 5/7/15.
//  Copyright (c) 2015 Jeffrey J. Early. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>
#import <GLNumericalModelingKit/GLOperationOptimizer.h>
#import "MomentTensorModels.h"
#import "MomentEllipseFileReader.h"
#import "EllipseErrorOperation.h"
#import "DrifterTracksFileReader.h"

int main(int argc, const char * argv[])
{
	
	@autoreleasepool {
		
		GLEquation *equation = [[GLEquation alloc] init];
		
		//		NSString *filePath =  @"/Users/jearly/Documents/LatMix/drifters/observations/griddedRho1Drifters.txt";
		//		DrifterTracksFileReader *trackReader = [[DrifterTracksFileReader alloc] initWithURL: [NSURL fileURLWithPath: filePath] equation: equation];
		//		MomentTensorModels *model = [[MomentTensorModels alloc] initWithXPositions: trackReader.x yPositions:trackReader.y time:trackReader.t];
		
		NSFileManager *fileManager = [[NSFileManager alloc] init];
		NSString *folderPath = @"/Users/jearly/Documents/LatMix/drifters/observations/griddedRho2DrifterDay6MomementEllipses/";
		//        NSString *folderPath = @"/Users/jearly/Documents/LatMix/drifters/synthetic/moment-ellipses/synthetic-diffusive/";
		NSArray *ellipseFiles = [fileManager contentsOfDirectoryAtPath: folderPath error: nil];
		
		NSMutableString *outputData = [NSMutableString stringWithFormat: @""];
#if ELLIPSE_ERROR_METHOD == 0
		NSString *outputFile = @"BestFitParameters_area_divergence_local_area.m";
#elif ELLIPSE_ERROR_METHOD == 1
		NSString *outputFile = @"BestFitParameters_area_divergence_total_area.m";
#endif
		NSUInteger i=1;
		
		for ( NSString *filename in ellipseFiles)
		{
			if ( ![filename.pathExtension isEqualToString: @"txt"]) continue;
			if ( [filename containsString: @"_BestFitSummary"]) continue;
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
			
			for (NSUInteger maxPoint=0; maxPoint < file.t.nDataPoints; maxPoint += 4)
			{
				GLFunction *a = [a]
				MomentTensorModels *models = [[MomentTensorModels alloc] initWithA: file.a b: file.b theta: file.angle time: file.t];
				NSArray *result = [models bestFitToDiffusivityModel];

				GLScalar *minError = result[0];
				GLScalar *minKappa = result[1];

				[outputData appendFormat: @"model1_error(%lu)=%g; model1_kappa(%lu)=%g; ", i, *(minError.pointerValue), i, *(minKappa.pointerValue)];
				NSLog(@"%@---diffusivity model total error: %f @ (kappa)=(%.4f)", filename.lastPathComponent, *(minError.pointerValue), *(minKappa.pointerValue));



				result = [models bestFitToStrainDiffusivityModel];

				minError = result[0];
				minKappa = result[1];
				GLScalar *minSigma = result[2];
				GLScalar *minTheta = result[3];

				[outputData appendFormat: @"model2_error(%lu)=%g; model2_kappa(%lu)=%g; model2_sigma(%lu)=%g; model2_theta(%lu)=%g;\n", i, *(minError.pointerValue), i, *(minKappa.pointerValue), i, *(minSigma.pointerValue), i, *(minTheta.pointerValue)];
				NSLog(@"%@---strain-diffusivity model total error: %f (kappa,sigma,theta)=(%.4f,%.3g,%.1f)", filename.lastPathComponent, *(minError.pointerValue), *(minKappa.pointerValue),*(minSigma.pointerValue),(*(minTheta.pointerValue))*180./M_PI);


				result = [models bestFitToVorticityStrainDiffusivityModelWithStartPoint:@[minKappa, minSigma, minTheta]];

				minError = result[0];
				minKappa = result[1];
				minSigma = result[2];
				minTheta = result[3];
				GLScalar *minZeta = result[4];

				[outputData appendFormat: @"model3_error(%lu)=%g; model3_kappa(%lu)=%g; model3_sigma(%lu)=%g; model3_theta(%lu)=%g; model3_zeta(%lu)=%g;\n", i, *(minError.pointerValue), i, *(minKappa.pointerValue), i, *(minSigma.pointerValue), i, *(minTheta.pointerValue), i, *(minZeta.pointerValue)];

				NSLog(@"%@---vorticity-strain-diffusivity model total error: %f (kappa,sigma,theta,zeta)=(%.4f,%.3g,%.1f,%.3g)", filename.lastPathComponent, *(minError.pointerValue), *(minKappa.pointerValue),*(minSigma.pointerValue),(*(minTheta.pointerValue))*180./M_PI, *(minZeta.pointerValue));

				i++;
			}
		}
		
		[outputData writeToFile: [NSString stringWithFormat: @"%@%@", folderPath, outputFile] atomically: YES encoding: NSUTF8StringEncoding error: nil];
	}
	return 0;
}
