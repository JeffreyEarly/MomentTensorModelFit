//
//  EllipseErrorOperation.h
//  MomentTensorModelFit
//
//  Created by Jeffrey Early on 11/15/13.
//  Copyright (c) 2013 Jeffrey J. Early. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>

#define ELLIPSE_ERROR_METHOD 0

@interface EllipseErrorOperation : GLVariableOperation

- (EllipseErrorOperation *) initWithParametersFromEllipseA: (NSArray *) paramsA ellipseB: (NSArray *) paramsB;

@end
