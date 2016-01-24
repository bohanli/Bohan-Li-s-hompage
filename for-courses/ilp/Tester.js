function test() {
    //  min c.x
    //  st  A.x = b
    //      l <= x <= u
    //  Some x are integer, and b >= 0
    var model = {};
    model.A = [[ 9, 7, 1, 0],
              [ 7, 20, 0, 1]];
    model.b =  [ 56, 70];
    model.c =  [-4,-9, 0, 0];
    model.m = 2;
    model.n = 4;
    // test.xLB = [0, 0, 0, 0];
    // test.xUB = [Infinity, Infinity, Infinity, Infinity];
    // test.xINT = [true, true, true, true];
    IlpSolver.solveILP(model);
    console.log(model.x, model.z);
}

console.log("Testing Branch and Bound");
test();