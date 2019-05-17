//https://github.com/processing/p5.js/blob/master/developer_docs/webgl_mode_architecture.md

let pg;
var program;
preload = function()
{
  program = loadShader('js/shaders/vert.shader', 'js/shaders/frag.shader');
};

function setup()
{
  createCanvas(400, 400, WEBGL);
  gl = canvas.getContext('webgl');
  pg = createGraphics(400, 400);
  pg.textSize(75);
  pg.background(0, 0);
  texture(pg);
  shader(program);
}

function draw()
{
  background(255);
  beginShape(TRIANGLE_STRIP);
  vertex(30, 75, 0);
  vertex(40, 20,100);
  vertex(50, 75,0);
  vertex(60, 20,0);
  vertex(70, 75,0);
  vertex(80, 20,0);
  vertex(200, 75,0);
  endShape();
}
