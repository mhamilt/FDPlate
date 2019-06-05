var screensize = (4 * $(window).width()) / 5;

let pg;
var shdr;
let zscl = 150000.0;
preload = function()
{
  shdr = loadShader('js/shaders/vert.shader', 'js/shaders/frag.shader');
};

function setup()
{
  var canvas = createCanvas(screensize, screensize, WEBGL);
  frameRate(30);
  // gl = canvas.getContext('webgl');
  canvas.parent('sketch-holder');
  pg = createGraphics(200, 200);
  pg.textSize(75);
  pg.background(0, 100);
  fill(255);
  texture(pg);
  shader(shdr);
  console.log(u[0]);
}


function draw()
{
  drawPlate();
}


function drawPlate()
{
  for (var i = 0; i < 10; ++i)
  {
    fd_update();
  }
  background(0);
  var scale = width / (Nx - 1);
  push();
  translate(0, 0, -200);
  rotateX(PI / 3);
  rotateZ(millis() * 0.0002);

  // rotateY(millis() /1000);
  translate(-width / 2, -height / 2, 0);
  for (var yi = 0; yi < Ny; ++yi)
  {
    beginShape(TRIANGLE_STRIP);
    for (var xi = 0; xi < Nx; ++xi)
    {
      var cp = (xi) + ((yi) * Nx); // current povar
      vertex(xi * scale, yi * scale, Math.floor(u[cp] * zscl));
      vertex(xi * scale, (yi + 1) * scale, Math.floor(u[cp + Ny] * zscl));
    }
    endShape()
  }

  pop();
}

function mousePressed()
{
  console.log(u[Nx / 2 + Ny / 2 * Nx] * zscl);

  var point = Math.floor(((Nx - 1) * mouseY / height) + ((Ny - 1) * Nx * mouseX / width))
  console.log(point);
  u1[point] += 0.0001;
}

function mouseDragged()
{
  var point = Math.floor(((Nx - 1) * mouseY / height) + ((Ny - 1) * Nx * mouseX / width))

  u1[point] += 0.00002;
}
