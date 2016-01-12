var IlpSolver = {
	// some constants 
	INFEASIBLE: 0, OPTIMAL: 1, UNBOUNDED: 2, MAXITERATIONS: Infinity,

	solveILP: function(ILPModel) {
		//  min cx
		//  s.t.  Ax = b
		//      x are non-negtive integers

		var maxNodes = Infinity, log = [];
		var currentBest = Infinity, currentBestX = new Array(ILPModel.n); // current OPT
		var mostFracIndex, mostFracValue, fracValue; // most fractional value

		
		/* branch & bound tree */
		// init the branch & bound tree
		var unsolvedLPs = new Array();
		ILPModel.solved = false;
		unsolvedLPs.push(ILPModel);
		var nodeCount = 0;

		// specify ILPModel
		var n = ILPModel.n;
		ILPModel.xLB = new Array(n);
		ILPModel.xUB = new Array(n);
		ILPModel.xINT = new Array(n); 
		for( i = 0; i < n; i++ ) ILPModel.xLB[i] = 0;
		for( i = 0; i < n; i++ ) ILPModel.xUB[i] = Infinity;
		for( i = 0; i < n; i++ ) ILPModel.xINT[i] = true;

		// iteratively process unsolved nodes on the branch & bound tree (BFS)
		while ( unsolvedLPs.length >= 1 ) {
			nodeCount += 1;
			model = unsolvedLPs.shift();
			
			// stop if nodeCount >= maxNodes
			if ( nodeCount >= maxNodes ) {
				unsolvedLPs = [];
				ILPModel.status = IlpSolver.INFEASIBLE;
				console.log("No solution!");
				return;
			}
			
			// solve the LP at this node
			IlpSolver.solveLP(model);
			if ( model.status == IlpSolver.INFEASIBLE )
				continue;
			
			// compare with current OPT
			if ( model.z > currentBest )
				continue;
			
			// check each value of the LP solution
			mostFracIndex = -1, mostFracValue = 0;
			for ( i = 0; i < model.n; i++ ) if (Math.abs(Math.floor(model.x[i]) - model.x[i]) > 0.0001) {
				// most fractional value 
				fracValue = Math.min( Math.abs(Math.floor(model.x[i]) - model.x[i]), Math.abs(Math.ceil (model.x[i]) - model.x[i]));
				if (fracValue > mostFracValue) { mostFracIndex = i, mostFracValue = fracValue; }
			}
			
			// check whether a fractional value is found
			if ( mostFracIndex == -1 ) {
				// if no fractional value found, terminate the branch and update current OPT
				if ( model.z < currentBest ) {
					currentBest = model.z;
					for ( i = 0; i < model.n; i++ )
						currentBestX[i] = model.x[i];
				}
			} else {
				// if some fractional value found, do branching again on this node
				// lower branch
				newBranch1 = IlpSolver.copyModel(model);
				newBranch1.xUB[mostFracIndex] = Math.floor(newBranch1.x[mostFracIndex])
				newBranch1.z = model.z;
				unsolvedLPs.push(newBranch1);
				
				// higher branch
				newBranch2 = IlpSolver.copyModel(model);
				newBranch2.xLB[mostFracIndex] = Math.ceil(newBranch2.x[mostFracIndex])
				newBranch2.z = model.z;
				unsolvedLPs.push(newBranch2);
			}
		}
	
		// check whether there is an OPT
		ILPModel.nodeCount = nodeCount;
		if (currentBest < Infinity) {
			ILPModel.x = currentBestX;
			ILPModel.z = currentBest;
			ILPModel.status = IlpSolver.OPTIMAL;
			console.log("Done!");
		} else {
			ILPModel.status = IlpSolver.INFEASIBLE;
			console.log("No solution!");
		}
	},


	solveLP : function(model) {
		// 	min c.x
		//  st  A.x = b
		//      l <= x <= u

		log = [];

		var BASIS = 0, NONBASIS_L = +1, NONBASIS_U = -1, ZERO = 0.000001;

		A = model.A; b = model.b; c = model.c;
		m = model.m; n = model.n;
		xLB = model.xLB; xUB = model.xUB;

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
			if ( iter >= IlpSolver.MAXITERATIONS ) {
				z = 0.0;
				for (i = 0; i < n; i++)
					z += c[i] * x[i];
				model.z = z;
				model.x = x;
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
			var minRC = -ZERO, s = -1;
			for ( i = 0; i < n; i++ ) if ( xStat[i] * rc[i] < minRC ) { 
				minRC = xStat[i] * rc[i];
				s = i; 
			}

			// if no entering variable
			if ( s == -1 ) {
				if ( isStageOne ) {
					z = 0.0;
					for ( i = 0; i < m; i++ )
						z += Cb[i] * x[basisVars[i]];

					if ( z > ZERO ) {
						model.status = IlpSolver.INFEASIBLE;
						break;
					} else { // state two
						isStageOne = false;
						for ( i = 0; i < m; i++ ) 
							Cb[i] = (basisVars[i] < n) ? (c[basisVars[i]]) : (0.0);
						continue;
					}
				} else {
					model.status = IlpSolver.OPTIMAL;
					z = 0.0;
					for (i = 0; i < n; i++)
						z += c[i] * x[i];
					model.z = z;
					model.x = x;
					break;
				}
			}
			
			// calculate BinvAs
			for ( i = 0; i < m; i++ ) {
				BinvAs[i] = 0.0;
				for ( k = 0; k < m; k++ )
					BinvAs[i] += Binv[i][k] * A[k][s];
			}

			// test ratio
			var minRatio = Infinity, ratio = 0.0, r = -1;
			var rIsEV = false;
			ratio = xUB[s] - xLB[s];
			if ( ratio <= minRatio ) {
				minRatio = ratio;
				r = -1;
				rIsEV = true;
			}

			for ( i = 0; i < m; i++ ) {
				j = basisVars[i];
				var jLB = (j >= n) ? 0.0 : xLB[j];
				var jUB = (j >= n) ? Infinity : xUB[j];

				if ( -1*xStat[s]*BinvAs[i] > +ZERO ) { 
					ratio = (x[j] - jUB) / (xStat[s]*BinvAs[i]);
					if ( ratio <= minRatio ) { 
						minRatio = ratio; 
						r = i; 
						rIsEV = false; 
					}
				}
				if ( +1*xStat[s]*BinvAs[i] > +ZERO ) {
					ratio = (x[j] - jLB) / (xStat[s]*BinvAs[i]);
					if ( ratio <= minRatio ) { 
						minRatio = ratio; 
						r = i; 
						rIsEV = false; 
					}
				}
			}
				
			if ( minRatio >= Infinity ) {
				if ( !isStageOne ) model.status = IlpSolver.UNBOUNDED;
				break;
			}
			
			// update solution and basis
			x[s] += xStat[s] * minRatio;
			for ( i = 0; i < m; i++ )
				x[basisVars[i]] -= xStat[s] * minRatio * BinvAs[i];

			if ( !rIsEV ) {
				var erBinvAs = BinvAs[r];

				for ( i = 0; i < m; i++ ) if (i != r) {
					var eiBinvAsOvererBinvAs = BinvAs[i] / erBinvAs;
					for ( j = 0; j < m; j++ ) 
						Binv[i][j] -= eiBinvAsOvererBinvAs * Binv[r][j]
				}

				for ( j = 0; j < m; j++ )
					Binv[r][j] /= erBinvAs;

				xStat[s] = BASIS;
				if ( basisVars[r] < n ) {
					if ( Math.abs(x[basisVars[r]] - xLB[basisVars[r]]) < ZERO ) xStat[basisVars[r]] = NONBASIS_L;
					if ( Math.abs(x[basisVars[r]] - xUB[basisVars[r]]) < ZERO ) xStat[basisVars[r]] = NONBASIS_U;
				} else {
					if ( Math.abs(x[basisVars[r]] - 0.00000) < ZERO ) xStat[basisVars[r]] = NONBASIS_L;
					if ( Math.abs(x[basisVars[r]] - Infinity) < ZERO ) xStat[basisVars[r]] = NONBASIS_U;
				}

				Cb[r] = isStageOne ? 0.0 : c[s];
				basisVars[r] = s;

			} else {
				// degenerate iteration
				if ( xStat[s] == NONBASIS_L ) xStat[s] = NONBASIS_U;
				else xStat[s] = NONBASIS_L;
			}
		}
	},


	copyModel: function(model) {
		// copyModel: deeply copy the model

		if (model == null || typeof(model) != 'object')
			return model;

		var newModel = new model.constructor();

		for (var key in model)
			newModel[key] = IlpSolver.copyModel(model[key]);

		return newModel;
	}

};