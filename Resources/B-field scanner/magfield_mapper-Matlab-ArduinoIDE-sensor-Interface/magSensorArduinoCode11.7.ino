const float xPin = A0;  // Analog pin for X-axis. Not sure if float is correct/will work. But we want the data outputed to be in decimals for fidelity
const float yPin = A1;  // Analog pin for Y-axis
const float zPin = A2;  // Analog pin for Z-axis
float initialXValue = 0.0;
float initialYValue = 0.0;
float initialZValue = 0.0;

char strX[20]; //at least as long as # of char +1
char strY[20]; //are these too long/short. Idk but seems to work/not matter
char strZ[20]; 
char buffer[50];

void setup() {
  Serial.begin(9600);
  pinMode(xPin, INPUT);
  pinMode(yPin, INPUT);
  pinMode(zPin, INPUT);
  initialXValue = analogRead(xPin);
  initialYValue = analogRead(yPin);
  initialZValue = analogRead(zPin);
}
void loop() {
  float xValue = analogRead(xPin) - initialXValue;
  float yValue = analogRead(yPin) - initialYValue;
  float zValue = analogRead(zPin) - initialZValue;
  dtostrf(xValue, -5, 7, strX); 
  dtostrf(yValue, -5, 7, strY);
  dtostrf(zValue, -5, 7, strZ);
  sprintf(buffer, "%s, %s, %s", strX, strY, strZ); //Creates single line string to output data. Useful for trying to format data when grabbing data in labview.
  Serial.println(buffer); //Data is outputed as string "XX.00000, YY.000000, ZZ.000000"
  delay(10);  // rate in millisenconds
}
