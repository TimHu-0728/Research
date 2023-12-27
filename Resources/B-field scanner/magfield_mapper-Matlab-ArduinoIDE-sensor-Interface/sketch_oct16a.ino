const int xPin = A0;  // Analog pin for X-axis
const int yPin = A1;  // Analog pin for Y-axis
const int zPin = A2;  // Analog pin for Z-axis
int initialXValue = 0;
int initialYValue = 0;
int initialZValue = 0;
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
  int xValue = analogRead(xPin) - initialXValue;
  int yValue = analogRead(yPin) - initialYValue;
  int zValue = analogRead(zPin) - initialZValue;
  Serial.print(xValue);
  Serial.println("");
  Serial.print(yValue);
  Serial.println("");
  Serial.print(zValue);
  Serial.println("");
  delay(1000);  // rate in millisenconds
}
