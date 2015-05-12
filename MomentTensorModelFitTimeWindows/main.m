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

// 0 is an extending time window
// 1 is a rollowing time window
#define WINDOWING 1

int main(int argc, const char * argv[])
{
	
	@autoreleasepool {
		
		GLEquation *equation = [[GLEquation alloc] init];
		
		//		NSString *filePath =  @"/Users/jearly/Documents/LatMix/drifters/observations/griddedRho1Drifters.txt";
		//		DrifterTracksFileReader *trackReader = [[DrifterTracksFileReader alloc] initWithURL: [NSURL fileURLWithPath: filePath] equation: equation];
		//		MomentTensorModels *model = [[MomentTensorModels alloc] initWithXPositions: trackReader.x yPositions:trackReader.y time:trackReader.t];
		
		NSFileManager *fileManager = [[NSFileManager alloc] init];
		NSString *folderPath = @"/Users/jearly/Documents/Models/LatMix/drifters/observations/griddedRho2DrifterMomementEllipses/";
		//        NSString *folderPath = @"/Users/jearly/Documents/LatMix/drifters/synthetic/moment-ellipses/synthetic-diffusive/";
		NSArray *ellipseFiles = [fileManager contentsOfDirectoryAtPath: folderPath error: nil];
		
		NSMutableString *outputData = [NSMutableString stringWithFormat: @""];
#if ELLIPSE_ERROR_METHOD == 0
	#if WINDOWING == 0
			[outputData appendFormat: @"titleText=\'Extending window in six hour increments, local area divergence\';"];
		NSString *outputFile = @"BestFitParams_area_div_local_area_extending_time_window.m";
	#elif WINDOWING == 1
		NSUInteger windowLength = 120;
        NSUInteger windowIncrement = 2;
		[outputData appendFormat: @"titleText=\'Rolling %lu hour window in %lu hour increments, local area divergence\';",windowLength/2,windowIncrement/2];
		NSString *outputFile = [NSString stringWithFormat:@"BestFitParams_area_div_local_area_rolling_%lu_hour_window.m",windowLength/2];
	#endif
#elif ELLIPSE_ERROR_METHOD == 1
	#if WINDOWING == 0
		[outputData appendFormat: @"titleText=\'Extending window in six hour increments\';"];
		NSString *outputFile = @"BestFitParams_area_div_total_area_extending_time_window.m";
	#elif WINDOWING == 1
		NSUInteger windowLength = 20;
        NSUInteger windowIncrement = 2;
		[outputData appendFormat: @"titleText=\'Rolling %lu hour window in %lu hour increments\';",windowLength/2,windowIncrement/2];
		NSString *outputFile = [NSString stringWithFormat:@"BestFitParams_area_div_total_area_rolling_%lu_hour_window.m",windowLength/2];
	#endif
#endif
		NSUInteger ensemble=1;
		
		for ( NSString *filename in ellipseFiles)
		{
			if ( ![filename.pathExtension isEqualToString: @"txt"]) continue;
			//if ( [filename containsString: @"_BestFitSummary"]) continue;
			if ( [[filename substringToIndex: 1] isEqualToString: @"_"]) continue;
			MomentEllipseFileReader *file = [[MomentEllipseFileReader alloc] initWithURL: [NSURL fileURLWithPath: [NSString stringWithFormat: @"%@%@", folderPath, filename]] equation: equation];
			if (!file) continue;
			
			NSArray *drifterIDs = [[filename stringByDeletingPathExtension] componentsSeparatedByString: @"_"];
			NSMutableString *idString = [NSMutableString stringWithFormat: @"drifterIDs{%lu}=[", ensemble];
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
			
			NSUInteger timeWindow = 1;
#if WINDOWING == 0
			
			NSUInteger minPoint = 0;
			for (NSUInteger windowLength=12; minPoint+windowLength < file.t.nDataPoints; windowLength += 12)
#elif WINDOWING == 1
			for (NSUInteger minPoint=0; minPoint+windowLength/2 < file.t.nDataPoints; minPoint += windowIncrement)
#endif
			{
				NSRange timeRange = NSMakeRange(minPoint, windowLength);
				if (timeRange.location + timeRange.length >= file.t.nDataPoints	) {
					timeRange.length = file.t.nDataPoints - 1 - timeRange.location;
				}
				NSValue *rangeValue = [NSValue valueWithRange: timeRange];
				GLFunction *a = [file.a variableFromIndexRange: @[rangeValue]];
				GLFunction *b = [file.b variableFromIndexRange: @[rangeValue]];
				GLFunction *angle = [file.angle variableFromIndexRange: @[rangeValue]];
				GLFunction *t = [file.t variableFromIndexRange: @[rangeValue]];
				GLFunction *t0 = [t minus:@(t.pointerValue[0])];
				MomentTensorModels *models = [[MomentTensorModels alloc] initWithA: a b: b theta: angle time: t0];
				
				NSArray *result = [models bestFitToDiffusivityModel];

				GLScalar *minError = result[0];
				GLScalar *minKappa = result[1];
				
				[outputData appendFormat: @"timeStart(%lu,%lu)=%g; timeEnd(%lu,%lu)=%g; ",timeWindow, ensemble, t.pointerValue[0],timeWindow, ensemble, t.pointerValue[timeRange.length-1]];
				NSLog(@"timeStart(%lu,%lu)=%g; timeEnd(%lu,%lu)=%g; ",timeWindow, ensemble, t.pointerValue[0],timeWindow, ensemble, t.pointerValue[timeRange.length-1]);
				
				[outputData appendFormat: @"model1_error(%lu,%lu)=%g; model1_kappa(%lu,%lu)=%g; ",timeWindow, ensemble, *(minError.pointerValue),timeWindow, ensemble, *(minKappa.pointerValue)];
				NSLog(@"%@---diffusivity model total error: %f @ (kappa)=(%.4f)", filename.lastPathComponent, *(minError.pointerValue), *(minKappa.pointerValue));



				result = [models bestFitToStrainDiffusivityModel];

				minError = result[0];
				minKappa = result[1];
				GLScalar *minSigma = result[2];
				GLScalar *minTheta = result[3];

				[outputData appendFormat: @"model2_error(%lu,%lu)=%g; model2_kappa(%lu,%lu)=%g; model2_sigma(%lu,%lu)=%g; model2_theta(%lu,%lu)=%g;\n",timeWindow, ensemble, *(minError.pointerValue),timeWindow, ensemble, *(minKappa.pointerValue),timeWindow, ensemble, *(minSigma.pointerValue),timeWindow, ensemble, *(minTheta.pointerValue)];
				NSLog(@"%@---strain-diffusivity model total error: %f (kappa,sigma,theta)=(%.4f,%.3g,%.1f)", filename.lastPathComponent, *(minError.pointerValue), *(minKappa.pointerValue),*(minSigma.pointerValue),(*(minTheta.pointerValue))*180./M_PI);


				result = [models bestFitToVorticityStrainDiffusivityModelWithStartPoint:@[minKappa, minSigma, minTheta]];

				minError = result[0];
				minKappa = result[1];
				minSigma = result[2];
				minTheta = result[3];
				GLScalar *minZeta = result[4];

				[outputData appendFormat: @"model3_error(%lu,%lu)=%g; model3_kappa(%lu,%lu)=%g; model3_sigma(%lu,%lu)=%g; model3_theta(%lu,%lu)=%g; model3_zeta(%lu,%lu)=%g;\n",timeWindow, ensemble, *(minError.pointerValue),timeWindow, ensemble, *(minKappa.pointerValue),timeWindow, ensemble, *(minSigma.pointerValue),timeWindow, ensemble, *(minTheta.pointerValue),timeWindow, ensemble, *(minZeta.pointerValue)];

				NSLog(@"%@---vorticity-strain-diffusivity model total error: %f (kappa,sigma,theta,zeta)=(%.4f,%.3g,%.1f,%.3g)", filename.lastPathComponent, *(minError.pointerValue), *(minKappa.pointerValue),*(minSigma.pointerValue),(*(minTheta.pointerValue))*180./M_PI, *(minZeta.pointerValue));
				timeWindow++;
			}
			ensemble++;
		}
		
		[outputData writeToFile: [NSString stringWithFormat: @"%@%@", folderPath, outputFile] atomically: YES encoding: NSUTF8StringEncoding error: nil];
	}
	return 0;
}
