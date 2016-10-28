//
//  main.m
//  MomentTensorModelFitFixedParameters
//
//  Created by  on 10/28/16.
//  Copyright © 2016 Jeffrey J. Early. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>
#import <GLNumericalModelingKit/GLOperationOptimizer.h>
#import "MomentTensorModels.h"
#import "MomentEllipseFileReader.h"
#import "EllipseErrorOperation.h"
#import "DrifterTracksFileReader.h"

// 0 is a single model fit to the whole time series.
// 1 is a rolling time window
#define WINDOWING 1
#define WINDOW_LENGTH_IN_HOURS 24
#define SITE 1

int main(int argc, const char * argv[])
{
    
    @autoreleasepool {
        
        GLEquation *equation = [[GLEquation alloc] init];
        
        NSFileManager *fileManager = [[NSFileManager alloc] init];
#if SITE == 1
        NSString *folderPath = @"/Users/jearly/Documents/ProjectRepositories/LatMix/drifters/observations/griddedRho1DrifterMomementEllipses/";
#elif SITE == 2
        NSString *folderPath = @"/Users/jearly/Documents/ProjectRepositories/LatMix/drifters/observations/griddedRho2DrifterMomementEllipses/";
#endif
        NSArray *ellipseFiles = [fileManager contentsOfDirectoryAtPath: folderPath error: nil];
        
        NSMutableString *outputData = [NSMutableString stringWithFormat: @""];
#if ELLIPSE_ERROR_METHOD == 0
#if WINDOWING == 0
        [outputData appendFormat: @"titleText=\'Full experiment fit\';"];
        NSString *outputFile = @"BestFitParams_area_div_local_area.m";
#elif WINDOWING == 1
        NSUInteger windowLength = WINDOW_LENGTH_IN_HOURS*4+1;
        NSUInteger windowIncrement = 4;
        [outputData appendFormat: @"titleText=\'Rolling %d hour window in %d hour increments, local area divergence\';",WINDOW_LENGTH_IN_HOURS,windowIncrement/4];
        NSString *outputFile = [NSString stringWithFormat:@"BestFitCustomParams_area_div_local_area_rolling_%d_hour_window.m",WINDOW_LENGTH_IN_HOURS];
#endif
#elif ELLIPSE_ERROR_METHOD == 1
#if WINDOWING == 0
        [outputData appendFormat: @"titleText=\'Full experiment fit\';"];
        NSString *outputFile = @"BestFitParams_area_div_total_area.m";
#elif WINDOWING == 1
        // New time series is in 15 minute increments, not 30.
        NSUInteger windowLength = WINDOW_LENGTH_IN_HOURS*4+1;
        NSUInteger windowIncrement = 4;
        [outputData appendFormat: @"titleText=\'Rolling %d hour window in %lu hour increments\';",WINDOW_LENGTH_IN_HOURS,windowIncrement/4];
        NSString *outputFile = [NSString stringWithFormat:@"BestFitCustomParams_area_div_total_area_rolling_%d_hour_window_FineGrid.m",WINDOW_LENGTH_IN_HOURS];
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
            NSUInteger windowLength = file.t.nDataPoints;
            for (NSUInteger minPoint=0; minPoint+windowLength <= file.t.nDataPoints; minPoint += 1)
#elif WINDOWING == 1
                for (NSUInteger minPoint=0; minPoint+windowLength <= file.t.nDataPoints; minPoint += windowIncrement)
#endif
                {
                    NSRange timeRange = NSMakeRange(minPoint, windowLength);
                    if (timeRange.location + timeRange.length >= file.t.nDataPoints	) {
                        timeRange.length = file.t.nDataPoints - timeRange.location;
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
                    
                    
                    result = [models bestFitToStrainDiffusivityModelWithFixedStrainAngle:[GLScalar scalarWithValue: -33.0*M_PI/180. forEquation: equation]];
                    
                    minError = result[0];
                    minKappa = result[1];
                    GLScalar *minSigma = result[2];
                    GLScalar *minTheta = result[3];
                    
                    [outputData appendFormat: @"model2_error(%lu,%lu)=%g; model2_kappa(%lu,%lu)=%g; model2_sigma(%lu,%lu)=%g; model2_theta(%lu,%lu)=%g;\n",timeWindow, ensemble, *(minError.pointerValue),timeWindow, ensemble, *(minKappa.pointerValue),timeWindow, ensemble, *(minSigma.pointerValue),timeWindow, ensemble, *(minTheta.pointerValue)];
                    NSLog(@"%@---strain-diffusivity model total error: %f (kappa,sigma,theta)=(%.4f,%.3g,%.1f)", filename.lastPathComponent, *(minError.pointerValue), *(minKappa.pointerValue),*(minSigma.pointerValue),(*(minTheta.pointerValue))*180./M_PI);
                    
                    
                    
                    
                    result = [models bestFitToVorticityStrainDominatedDiffusivityModelWithStartPoint:@[minKappa, minSigma, minTheta]];
                    
                    minError = result[0];
                    minKappa = result[1];
                    minSigma = result[2];
                    minTheta = result[3];
                    GLScalar *minZeta = result[4];
                    
                    [outputData appendFormat: @"model3_error(%lu,%lu)=%g; model3_kappa(%lu,%lu)=%g; model3_sigma(%lu,%lu)=%g; model3_theta(%lu,%lu)=%g; model3_zeta(%lu,%lu)=%g;\n",timeWindow, ensemble, *(minError.pointerValue),timeWindow, ensemble, *(minKappa.pointerValue),timeWindow, ensemble, *(minSigma.pointerValue),timeWindow, ensemble, *(minTheta.pointerValue),timeWindow, ensemble, *(minZeta.pointerValue)];
                    
                    NSLog(@"%@---vorticity-strain(dominate)-diffusivity model total error: %f (kappa,sigma,theta,zeta)=(%.4f,%.3g,%.1f,%.3g)", filename.lastPathComponent, *(minError.pointerValue), *(minKappa.pointerValue),*(minSigma.pointerValue),(*(minTheta.pointerValue))*180./M_PI, *(minZeta.pointerValue));
                    
                    
                    
                    
                    result = [models bestFitToVorticityDiffusivityModel];
                    
                    minError = result[0];
                    minKappa = result[1];
                    minZeta = result[2];
                    
                    [outputData appendFormat: @"model4_error(%lu,%lu)=%g; model4_kappa(%lu,%lu)=%g; model4_zeta(%lu,%lu)=%g;\n",timeWindow, ensemble, *(minError.pointerValue),timeWindow, ensemble, *(minKappa.pointerValue),timeWindow, ensemble, *(minZeta.pointerValue)];
                    NSLog(@"%@---vorticity-diffusivity model total error: %f (kappa,zeta)=(%.4f,%.3g)", filename.lastPathComponent, *(minError.pointerValue), *(minKappa.pointerValue),*(minZeta.pointerValue));
                    
                    
                    
                    result = [models bestFitToVorticityStrainMatchedDiffusivityModelWithStartPoint:@[minKappa, minSigma, minTheta]];
                    
                    minError = result[0];
                    minKappa = result[1];
                    minSigma = result[2];
                    minTheta = result[3];
                    minZeta = result[4];
                    
                    [outputData appendFormat: @"model5_error(%lu,%lu)=%g; model5_kappa(%lu,%lu)=%g; model5_sigma(%lu,%lu)=%g; model5_theta(%lu,%lu)=%g; model5_zeta(%lu,%lu)=%g;\n",timeWindow, ensemble, *(minError.pointerValue),timeWindow, ensemble, *(minKappa.pointerValue),timeWindow, ensemble, *(minSigma.pointerValue),timeWindow, ensemble, *(minTheta.pointerValue),timeWindow, ensemble, *(minZeta.pointerValue)];
                    
                    NSLog(@"%@---vorticity-strain-diffusivity-matched model total error: %f (kappa,sigma,theta,zeta)=(%.4f,%.3g,%.1f,%.3g)", filename.lastPathComponent, *(minError.pointerValue), *(minKappa.pointerValue),*(minSigma.pointerValue),(*(minTheta.pointerValue))*180./M_PI, *(minZeta.pointerValue));
                    
                    
                    
                    result = [models bestFitToVorticityDominatedStrainDiffusivityModelWithStartPoint:@[minKappa, minSigma, minTheta]];
                    
                    minError = result[0];
                    minKappa = result[1];
                    minSigma = result[2];
                    minTheta = result[3];
                    minZeta = result[4];
                    
                    [outputData appendFormat: @"model6_error(%lu,%lu)=%g; model6_kappa(%lu,%lu)=%g; model6_sigma(%lu,%lu)=%g; model6_theta(%lu,%lu)=%g; model6_zeta(%lu,%lu)=%g;\n",timeWindow, ensemble, *(minError.pointerValue),timeWindow, ensemble, *(minKappa.pointerValue),timeWindow, ensemble, *(minSigma.pointerValue),timeWindow, ensemble, *(minTheta.pointerValue),timeWindow, ensemble, *(minZeta.pointerValue)];
                    
                    NSLog(@"%@---vorticity(dominant)-strain-diffusivity model total error: %f (kappa,sigma,theta,zeta)=(%.4f,%.3g,%.1f,%.3g)", filename.lastPathComponent, *(minError.pointerValue), *(minKappa.pointerValue),*(minSigma.pointerValue),(*(minTheta.pointerValue))*180./M_PI, *(minZeta.pointerValue));
                    
                    
                    timeWindow++;
                }
            ensemble++;
        }
        
        [outputData writeToFile: [NSString stringWithFormat: @"%@%@", folderPath, outputFile] atomically: YES encoding: NSUTF8StringEncoding error: nil];
    }
    return 0;
}

