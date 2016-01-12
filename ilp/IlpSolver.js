var IlpSolver = {
	// some constants 
	INFEASIBLE: 0, OPTIMAL: 1, UNBOUNDED: 2, MAXITERATIONS: 100,  ZERO: 0.000001, INF: 10000000,

	solveILP: function(ILPModel) {
		//  min cx
		//  s.t.  Ax = b
		//      x are non-negtive integers

		var currentBest = this.INF, currentBestX = new Array(ILPModel.n); // current OPT
		var mostFracIndex, mostFracValue, fracValue; // most fractional value
		var log = [];
		
		/* branch & bound tree */
		// init the branch & bound tree
		var unsolvedLPs = new Array();
		ILPModel.solved = false;
		unsolvedLPs.push(ILPModel);
		var nodeCnt = 0;

		// specify ILPModel
		var n = ILPModel.n;
		ILPModel.xLB = new Array(n);
		ILPModel.xUB = new Array(n);
		//ILPModel.xINT = new Array(n); 
		for( i = 0; i < n; i++ ) ILPModel.xLB[i] = 0;
		for( i = 0; i < n; i++ ) ILPModel.xUB[i] = this.INF;
		//for( i = 0; i < n; i++ ) ILPModel.xINT[i] = true;
		console.log(ILPModel);
		console.log("test");

		// iteratively process unsolved nodes on the branch & bound tree (BFS)
		while( unsolvedLPs.length >= 1 ) {
			nodeCnt += 1;
			model = unsolvedLPs.shift();
			
			// stop if nodeCnt >= MAXITERATIONS
			if( nodeCnt >= this.MAXITERATIONS ) {
				unsolvedLPs = [];
				ILPModel.status = this.INFEASIBLE;
				console.log("No solution!");
				return;
			}
			
			// solve the LP at this node
			console.log("=========================")
			console.log("unsolvedLPs", unsolvedLPs.length)
			console.log("model.old.xLB", model.xLB)
			console.log("model.old.xUB", model.xUB)
			this.solveLP(model);
			if( model.status == this.INFEASIBLE ) {
				continue;
			}
				
			console.log(model.x)

			// compare with current OPT
			console.log("model.z", model.z)
			console.log("currentBest", currentBest)
			if( model.z > currentBest )
				continue;

			console.log("model.z", model.z)
			
			// check each value of the LP solution

			
			console.log("model", model)
			mostFracIndex = -1, mostFracValue = 0;
			for( i = 0; i < model.n; i++ ) {
				if( model.x[i] < model.xLB[i] || model.x[i] > model.xUB[i] ) {
					continue;
				}


				console.log("model.x["+i+"]", model.x[i])
				console.log("Math.abs(Math.floor(model.x[i]) - model.x[i])", Math.abs(Math.floor(model.x[i]) - model.x[i]))
				if( Math.abs(Math.floor(model.x[i]) - model.x[i]) > this.ZERO ) {
					// most fractional value 
					//console.log("fuck")

					fracValue = Math.min( Math.abs(Math.floor(model.x[i]) - model.x[i]), Math.abs(Math.ceil (model.x[i]) - model.x[i]));
					if( fracValue > mostFracValue ) {
						//console.log("fracValue", fracValue)
						mostFracIndex = i;
						mostFracValue = fracValue; 
					}
				}

			}
			
			console.log("mostFracIndex", mostFracIndex)

			// check whether a fractional value is found
			if( mostFracIndex == -1 ) {
				// if no fractional value found, terminate the branch and update current OPT
				console.log("if( mostFracIndex == -1 ) ", model.z);
				if( model.z < currentBest && Math.abs(Math.floor(model.z) - model.z) < this.ZERO) {
					console.log("currentBest", currentBest)
					console.log("currentBestX", currentBestX)
					currentBest = model.z;
					for( i = 0; i < model.n; i++ )
						currentBestX[i] = model.x[i];
					console.log("NEW currentBest", currentBest)
					console.log("NEW currentBestX", currentBestX)
				}
			} else if( model.z < currentBest ) {

				console.log("model xLB[mostFracIndex]", model.xLB[mostFracIndex])
				console.log("model xUB[mostFracIndex]", model.xUB[mostFracIndex])

				// if some fractional value found, do branching again on this node
				// lower branch
				newBranch1 = this.copyModel(model);
				newBranch1.xUB[mostFracIndex] = Math.floor(model.x[mostFracIndex])


				console.log("lower branch xLB[mostFracIndex]", newBranch1.xLB[mostFracIndex])
				console.log("lower branch xUB[mostFracIndex]", newBranch1.xUB[mostFracIndex])
				// newBranch1.z = model.z;
				if( newBranch1.xLB[mostFracIndex] <= newBranch1.xUB[mostFracIndex]  )
					unsolvedLPs.push(newBranch1);
				
				// upper branch
				newBranch2 = this.copyModel(model);
				newBranch2.xLB[mostFracIndex] = Math.ceil(model.x[mostFracIndex])

				// if( newBranch2.xLB[mostFracIndex] > newBranch2.xUB[mostFracIndex] )
				// 	break;

				console.log("upper branch xLB[mostFracIndex]", newBranch2.xLB[mostFracIndex])
				console.log("upper branch xUB[mostFracIndex]", newBranch2.xUB[mostFracIndex])
				// newBranch2.z = model.z;
				if( newBranch2.xLB[mostFracIndex] <= newBranch2.xUB[mostFracIndex]  )
					unsolvedLPs.push(newBranch2);
			}
		}
	
		// check whether there is an OPT
		this.nodeCnt = nodeCnt;
		if( currentBest < this.INF ) {
			ILPModel.x = currentBestX;
			ILPModel.z = currentBest;
			ILPModel.status = this.OPTIMAL;
			console.log("Done!");
		} else {
			ILPModel.status = this.INFEASIBLE;
			console.log("No solution!");
		}

		for( i = 0; i < ILPModel.n; i++ ) 
				ILPModel.x[i] = Math.floor(ILPModel.x[i]);

		// for( i = 0; i < ILPModel.n; i++ ) {
		// 	if( Math.abs(Math.floor(ILPModel.x[i]) - ILPModel.x[i]) > this.ZERO ) {
		// 		ILPModel.x[i] = Math.floor(ILPModel.x[i]);
		// 	}

		// }

	},


	solveLP : function(model) {
		// 	min c.x
		//  st  A.x = b
		//      l <= x <= u

		var BASIS = 0, NONBASIS_L = +1, NONBASIS_U = -1;
		A = model.A; b = model.b; c = model.c; m = model.m; n = model.n; xLB = model.xLB; xUB = model.xUB;

		var xStat = new Array(n + m);
		var basisVars = new Array(m);
		var Binv = new Array(m); 
		var Cb = new Array(m);
		var P = new Array(m);
		var rc = new Array(n);
		var BinvAs = new Array(m);

		for ( i = 0; i < m; i++ )
			Binv[i] = new Array(m);

		var x = new Array(n + m), z, status; // solution
		
		/* init the initial solution */
		// original
		for ( i = 0; i < n; i++ ) {
			var absLB = Math.abs(xLB[i]), absUB = Math.abs(xUB[i]);
			x[i] = (absLB < absUB) ? xLB[i] : xUB[i];
			xStat[i] = (absLB < absUB) ? NONBASIS_L : NONBASIS_U;
		}

		// artificial (initial basis)
		for ( i = 0; i < m; i++ ) {
			x[i+n] = b[i];
			for ( j = 0; j < n; j++ )
				x[i+n] -= A[i][j] * x[j];

			xStat[i+n] = BASIS;
			basisVars[i] = i+n;

			Cb[i] = +1.0;
			for ( j = 0; j < m; j++ )
				Binv[i][j] = (i == j) ? 1.0 : 0.0;
		}

		
		/* iteratively solving Simplex */
		var isStageOne = true, iter = 0;
		while ( true ) {
			iter++;
			if ( iter >= this.MAXITERATIONS ) {
				z = 0.0;
				for (i = 0; i < n; i++)
					z += c[i] * x[i];
				model.z = z;
				model.x = x;
				//console.log("break")
				break;
			}

			for (i = 0; i < m; i++) {
				P[i] = 0.0;
				for (j = 0; j < m; j++) {
					P[i] += Cb[j] * Binv[j][i]
				}
			}

			for (j = 0; j < n; j++) {
				rc[j] = isStageOne ? 0.0 : c[j];
				for (i = 0; i < m; i++) {
					rc[j] -= P[i] * A[i][j];
				}
			}

			// pick entering variable
			var minRC = -this.ZERO, s = -1;
			for ( i = 0; i < n; i++ ) if ( xStat[i] * rc[i] < minRC ) { 
				minRC = xStat[i] * rc[i];
				s = i; 
			}

			//console.log("minRC", minRC)
			// if no entering variable
			if( s == -1 ) {
				if( isStageOne ) {
					z = 0.0;
					for ( i = 0; i < m; i++ )
						z += Cb[i] * x[basisVars[i]];

					if ( z > this.ZERO ) {
						model.status = this.INFEASIBLE;
						//console.log("break")
						break;
					} else { // state two
						isStageOne = false;
						for ( i = 0; i < m; i++ ) 
							Cb[i] = (basisVars[i] < n) ? (c[basisVars[i]]) : (0.0);
						continue;
					}
				} else {
					model.status = this.OPTIMAL;
					z = 0.0;
					for(i = 0; i < n; i++)
						z += c[i] * x[i];
					model.z = z;
					model.x = x;
					//console.log("break")
					break;
				}
			}
			
			// calculate BinvAs
			for( i = 0; i < m; i++ ) {
				BinvAs[i] = 0.0;
				for ( k = 0; k < m; k++ )
					BinvAs[i] += Binv[i][k] * A[k][s];
			}

			// test ratio
			var minRatio = Infinity, ratio = 0.0, r = -1;
			var rIsEV = false;
			ratio = xUB[s] - xLB[s];
			if( ratio <= minRatio ) {
				minRatio = ratio;
				r = -1;
				rIsEV = true;
			}


			// for( i = 0; i < m; i++ ) {

			// 	if( -1*xStat[s]*BinvAs[i] > +this.ZERO ) { 
				
			// 	}
			// }

			for( i = 0; i < m; i++ ) {
				j = basisVars[i];
				var jLB = (j >= n) ? 0.0 : xLB[j];
				var jUB = (j >= n) ? Infinity : xUB[j];

				if( -1*xStat[s]*BinvAs[i] > +this.ZERO ) { 
					ratio = (x[j] - jUB) / (xStat[s]*BinvAs[i]);
					if ( ratio <= minRatio ) { 
						minRatio = ratio; 
						r = i; 
						rIsEV = false; 
					}
				}
				if( +1*xStat[s]*BinvAs[i] > +this.ZERO ) {
					ratio = (x[j] - jLB) / (xStat[s]*BinvAs[i]);
					if ( ratio <= minRatio ) { 
						minRatio = ratio; 
						r = i; 
						rIsEV = false; 
					}
				}
			}
				
			if( minRatio >= Infinity ) {
				if ( !isStageOne ) model.status = this.UNBOUNDED;
				//console.log("break")
				break;
			}
			
			// update solution and basis
			x[s] += xStat[s] * minRatio;
			for( i = 0; i < m; i++ )
				x[basisVars[i]] -= xStat[s] * minRatio * BinvAs[i];

			if( !rIsEV ) {
				var erBinvAs = BinvAs[r];

				for( i = 0; i < m; i++ ) if (i != r) {
					var eiBinvAsOvererBinvAs = BinvAs[i] / erBinvAs;
					for ( j = 0; j < m; j++ ) 
						Binv[i][j] -= eiBinvAsOvererBinvAs * Binv[r][j]
				}

				for( j = 0; j < m; j++ )
					Binv[r][j] /= erBinvAs;

				xStat[s] = BASIS;
				if( basisVars[r] < n ) {
					if ( Math.abs(x[basisVars[r]] - xLB[basisVars[r]]) < this.ZERO ) xStat[basisVars[r]] = NONBASIS_L;
					if ( Math.abs(x[basisVars[r]] - xUB[basisVars[r]]) < this.ZERO ) xStat[basisVars[r]] = NONBASIS_U;
				} else {
					if ( Math.abs(x[basisVars[r]] - 0.00000) < this.ZERO ) xStat[basisVars[r]] = NONBASIS_L;
					if ( Math.abs(x[basisVars[r]] - Infinity) < this.ZERO ) xStat[basisVars[r]] = NONBASIS_U;
				}

				Cb[r] = isStageOne ? 0.0 : c[s];
				basisVars[r] = s;

			} else {
				// degenerate iteration
				if( xStat[s] == NONBASIS_L ) xStat[s] = NONBASIS_U;
				else xStat[s] = NONBASIS_L;
			}
		}
	},


	copyModel: function(model) {
		// copyModel: deeply copy the model

		if(model == null || typeof(model) != 'object')
			return model;

		var newModel = new model.constructor();

		for(var key in model)
			newModel[key] = this.copyModel(model[key]);

		return newModel;
	}

};