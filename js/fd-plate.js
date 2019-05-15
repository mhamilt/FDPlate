// Javascript FD Plate
//==============================================================================
//==============================================================================
// Constants
let pi = 3.14159265358979323846;
//==============================================================================
// Quick Signum Function used in the initial force
function signum(d)
{
  return (d <= 0) ? 0 : 1;
}
//==============================================================================
//----------------------------------------------------------------------------
// START EDIT HERE
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Flags
//----------------------------------------------------------------------------
// Conditions
var bctype = false; // boundary condition type: 0: simply supported, 1: clamped
var outtype = true; // output type: 0: displacement, 1: velocity

//----------------------------------------------------------------------------
// Parameters
//----------------------------------------------------------------------------

// simulation
var Tf = 2; // duration
var nu = .5; // Poisson Ratios (< .5)
var ctr[2] = {
  .45,
  .45
}; // centre povar of excitation as percentage
var wid = .25; // width (m)
var u0 = 0;
var v0 = 1; // excitation displacement and velocity
var rp[4] = {
  .45,
  .65,
  .85,
  .15
}; // readout position as percentage on grid.

//----------------------------------------------------------------------------
// Physical parameters
// // wood
var E = 11e9; // Young's modulus
var rho = 480; // density (kg/m^3)

// // steel
//	E = 2e11;							// Young's modulus
//	rho = 7850;							// density (kg/m^3)

var H = .005; // thickness (m)
var Lx = 1; // x-axis plate length (m)
var Ly = 1; // y-axis plate length (m)
var loss[4] = {
  100,
  8,
  1000,
  1
}; // loss [freq.(Hz), T60;...]

// I/O
var OSR = 1; // Oversampling ratio
var SR = 44.1e3; // sample rate (Hz)

//----------------------------------------------------------------------------
// END EDIT HERE
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Derived Parameters
//----------------------------------------------------------------------------
// Redefine SR from OSR
SR *= OSR;

// Motion Coefficients
var D = (E * (Math.pow(H, 3))) / (12 * (1 - Math.pow(nu, 2)));
var kappa = Math.sqrt(D / (rho * H));

//----------------------------------------------------------------------------
// Loss coefficients
//----------------------------------------------------------------------------
var z1 = 2 * kappa * (2 * pi * loss[0]) / (2 * Math.pow(kappa, 2));
var z2 = 2 * kappa * (2 * pi * loss[2]) / (2 * Math.pow(kappa, 2));
var sigma0 = 6 * Math.log(10) * (-z2 / loss[1] + z1 / loss[3]) / (z1 - z2);
var sigma1 = 6 * Math.log(10) * (1 / loss[1] - 1 / loss[3]) / (z1 - z2);

//----------------------------------------------------------------------------
// Grid Spacing
//----------------------------------------------------------------------------

var k = 1 / SR; // time step

// stability condition
var hmin = (Math.sqrt(4 * k * (sigma1 + Math.sqrt(Math.pow(sigma1, 2) + Math.pow(kappa, 2)))));

var Nx = Math.floor(Lx / hmin); // number of segments x-axis
var Ny = Math.floor(Ly / hmin); // number of segments y-axis
var h = Math.sqrt(Lx * Ly / (Nx * Ny));; // adjusted grid spacing x/y
Nx = Nx + 1;
Ny = Ny + 1; // grid povar number x and y
var mu = (kappa * k) / Math.pow(h, 2); // scheme parameter
var Nf = Math.floor(SR * Tf); // number of time steps
var ss = Nx * Ny; // total grid size.

//----------------------------------------------------------------------------
// Allocate Memory
//----------------------------------------------------------------------------
var u = new Array(ss);
var u1 = new Array(ss);
var u2 = new Array(ss);
var out = new Array(Nf);

//----------------------------------------------------------------------------
// Read In/Out
//----------------------------------------------------------------------------
var li = (Ny * (((ctr[1] * Nx) - 1))) + (ctr[0] * Ny) - 1;
var lo = (Ny * (((rp[1] * Nx) - 1))) + (rp[0] * Ny) - 1;

//----------------------------------------------------------------------------
// Excitation Force
//----------------------------------------------------------------------------
// raised Math.cosine in 2D
for (var xi = 1; xi < Nx - 1; xi++)
{
  var X = xi * h;

  for (var yi = 1; yi < Ny - 1; yi++)
  {
    var cp = yi + (xi * Ny);
    var Y = yi * h;
    var dist = Math.sqrt(Math.pow(X - (ctr[0] * Lx), 2) + Math.pow(Y - (ctr[1] * Ly), 2));
    var ind = signum((wid * 0.5) - dist); // displacement (Math.logical)
    var rc = .5 * ind * (1 + Math.cos(2 * pi * dist / wid)); // displacement
    u2[cp] = u0 * rc;
    u1[cp] = v0 * k * rc;
  }
}

//----------------------------------------------------------------------------
// Coefficient Matrices
//----------------------------------------------------------------------------

// coefficients are named based on position on the x and y axes.
var A00 = 1 / (1 + k * sigma0); // Central Loss Coeeffient (INVERTED)

//// Current time step (B) coeffients
// There are six unique coefficients for B coefs
var B00 = (-Math.pow(mu, 2) * 20 + (2 * sigma1 * k / Math.pow(h, 2)) * -4 + 2) * A00; // center
var B01 = (-Math.pow(mu, 2) * -8 + (2 * sigma1 * k / Math.pow(h, 2))) * A00; // 1-off
var B11 = (-Math.pow(mu, 2) * 2) * A00; // diag
var B02 = (-Math.pow(mu, 2) * 1) * A00; // 2-off

// Boundary Coefficients
var BC1, BC2;

if (bctype)
{
  BC1 = (-Math.pow(mu, 2) * 21 + (2 * sigma1 * k / Math.pow(h, 2)) * -4 + 2) * A00; // Side
  BC2 = (-Math.pow(mu, 2) * 22 + (2 * sigma1 * k / Math.pow(h, 2)) * -4 + 2) * A00; // Corner
}
else
{
  BC1 = (-Math.pow(mu, 2) * 19 + (2 * sigma1 * k / Math.pow(h, 2)) * -4 + 2) * A00; // Side
  BC2 = (-Math.pow(mu, 2) * 18 + (2 * sigma1 * k / Math.pow(h, 2)) * -4 + 2) * A00; // Corner
}

// Previous time step (C) coeffients
var C00 = (-(2 * sigma1 * k / Math.pow(h, 2)) * -4 - (1 - sigma0 * k)) * A00;
var C01 = -(2 * sigma1 * k / Math.pow(h, 2)) * A00;

//----------------------------------------------------------------------------
// Print Scheme Info
//----------------------------------------------------------------------------
console.Math.log("--- Coefficient Info --- \n\n");
console.Math.log("Loss A		: %.4fm \n", A00);
console.Math.log("Centre B    : %.4fm \n", B00);
console.Math.log("1-Grid B    : %.4fm \n", B01);
console.Math.log("2-Grid B	: %.4fm \n", B02);
console.Math.log("Diagonal B  : %.4fm \n", B11);
console.Math.log("Centre C	: %.4fm \n", C00);
console.Math.log("1-Grid C    : %.4fm \n", C01);
console.Math.log("Side Bound	: %.4fm \n", BC1);
console.Math.log("Cornr Bound : %.4fm \n", BC2);

console.Math.log("\n--- Scheme Info --- \n\n");
console.Math.log("Size		: %.1fm2 \n", Nx * h * Ny * h);
console.Math.log("Grid X-Ax   : %d \n", Nx);
console.Math.log("Grid Y-Ax   : %d \n", Ny);
console.Math.log("Total P		: %d \n", ss);
console.Math.log("Dur(samps)	: %d \n", Nf);
console.Math.log("In_cell		: %d\n", li);
console.Math.log("Out_cell	: %d\n", lo);
console.Math.log("Youngs		: %.2e\n", E);
console.Math.log("Sigma 0		: %f\n", sigma0);
console.Math.log("Sigma 1		: %f\n", sigma1);

//----------------------------------------------------------------------------
// Main Loop
//----------------------------------------------------------------------------

for (var n = 0; n < Nf; n++)
{
  for (var xi = 2; xi < Nx - 2; ++xi)
  {
    for (var yi = 2; yi < Ny - 2; ++yi)
    {
      var cp = (yi) + ((xi) * Ny); // current povar

      u[cp] = B00 * u1[cp] +
        B01 * (u1[cp - 1] + u1[cp + 1] + u1[cp - Ny] + u1[cp + Ny]) +
        B02 * (u1[cp - 2] + u1[cp + 2] + u1[cp - (2 * Ny)] + u1[cp + (2 * Ny)]) +
        B11 * (u1[cp - 1 - Ny] + u1[cp + 1 - Ny] + u1[cp + 1 + Ny] + u1[cp - 1 + Ny]) +
        C00 * u2[cp] +
        C01 * (u2[cp - 1] + u2[cp + 1] + u2[cp - Ny] + u2[cp + Ny]);
    }
  }

  // Update Side Boundaries
  //X-Axis

  for (var xi = 2; xi < Nx - 2; ++xi)
  {
    //North

    var cp = 1 + (xi * Ny); // current povar
    u[cp] = BC1 * u1[cp] +
      B01 * (u1[cp + 1] + u1[cp - Ny] + u1[cp + Ny]) +
      B02 * (u1[cp - 2] + u1[cp + 2] + u1[cp - (2 * Ny)] + u1[cp + (2 * Ny)]) +
      B11 * (u1[cp + 1 - Ny] + u1[cp + 1 + Ny]) +
      C00 * u2[cp] +
      C01 * (u2[cp + 1] + u2[cp - Ny] + u2[cp + Ny]);

    //South

    cp = Ny - 2 + (xi * Ny); // current povar
    u[cp] = BC1 * u1[cp] +
      B01 * (u1[cp - 1] + u1[cp - Ny] + u1[cp + Ny]) +
      B02 * (u1[cp - 2] + u1[cp - (2 * Ny)] + u1[cp + (2 * Ny)]) +
      B11 * (u1[cp - 1 - Ny] + u1[cp - 1 + Ny]) +
      C00 * u2[cp] +
      C01 * (u2[cp - 1] + u2[cp - Ny] + u2[cp + Ny]);


  }

  // Y-Axis

  for (var yi = 2; yi < Ny - 2; ++yi)
  {
    //West

    var cp = yi + Ny; // current povar
    u[cp] = BC1 * u1[cp] +
      B01 * (u1[cp - 1] + u1[cp + 1] + u1[cp + Ny]) +
      B02 * (u1[cp - 2] + u1[cp + 2] + u1[cp + (2 * Ny)]) +
      B11 * (u1[cp + 1 + Ny] + u1[cp - 1 + Ny]) +
      C00 * u2[cp] +
      C01 * (u2[cp - 1] + u2[cp + 1] + u2[cp + Ny]);

    //East

    cp = (yi) + Ny * (Nx - 2); // current povar
    u[cp] = BC1 * u1[cp] +
      B01 * (u1[cp - 1] + u1[cp + 1] + u1[cp - Ny]) +
      B02 * (u1[cp - 2] + u1[cp + 2] + u1[cp - (2 * Ny)]) +
      B11 * (u1[cp - 1 - Ny] + u1[cp + 1 - Ny]) +
      C00 * u2[cp] +
      C01 * (u2[cp - 1] + u2[cp + 1] + u2[cp - Ny]);

  }

  // Corner Boundaries
  {
    var cp = Ny + 1;
    u[cp] = BC2 * u1[cp] +
      B01 * (u1[cp - 1] + u1[cp + 1] + u1[cp - Ny] + u1[cp + Ny]) +
      B02 * (u1[cp + 2] + u1[cp + (2 * Ny)]) +
      B11 * (u1[cp - 1 - Ny] + u1[cp + 1 - Ny] + u1[cp + 1 + Ny] + u1[cp - 1 + Ny]) +
      C00 * u2[cp] +
      C01 * (u2[cp - 1] + u2[cp + 1] + u2[cp - Ny] + u2[cp + Ny]);

    cp = 2 * (Ny - 1);
    u[cp] = BC2 * u1[cp] +
      B01 * (u1[cp - 1] + u1[cp + 1] + u1[cp - Ny] + u1[cp + Ny]) +
      B02 * (u1[cp - 2] + u1[cp + (2 * Ny)]) +
      B11 * (u1[cp - 1 - Ny] + u1[cp + 1 - Ny] + u1[cp + 1 + Ny] + u1[cp - 1 + Ny]) +
      C00 * u2[cp] +
      C01 * (u2[cp - 1] + u2[cp + 1] + u2[cp - Ny] + u2[cp + Ny]);

    cp = Ny * (Nx - 2) + 1;
    u[cp] = BC2 * u1[cp] +
      B01 * (u1[cp - 1] + u1[cp + 1] + u1[cp - Ny] + u1[cp + Ny]) +
      B02 * (u1[cp + 2] + u1[cp - (2 * Ny)]) +
      B11 * (u1[cp - 1 - Ny] + u1[cp + 1 - Ny] + u1[cp + 1 + Ny] + u1[cp - 1 + Ny]) +
      C00 * u2[cp] +
      C01 * (u2[cp - 1] + u2[cp + 1] + u2[cp - Ny] + u2[cp + Ny]);

    cp = Ny * (Nx - 1) - 2;
    u[cp] = BC2 * u1[cp] +
      B01 * (u1[cp - 1] + u1[cp + 1] + u1[cp - Ny] + u1[cp + Ny]) +
      B02 * (u1[cp - 2] + u1[cp - (2 * Ny)]) +
      B11 * (u1[cp - 1 - Ny] + u1[cp + 1 - Ny] + u1[cp + 1 + Ny] + u1[cp - 1 + Ny]) +
      C00 * u2[cp] +
      C01 * (u2[cp - 1] + u2[cp + 1] + u2[cp - Ny] + u2[cp + Ny]);
  }

  out[n] = (outtype) ? SR * (u[lo] - u1[lo]) : u[lo];

  var dummy_ptr = u2;
  u2 = u1;
  u1 = u;
  u = dummy_ptr; // swap povarers
}