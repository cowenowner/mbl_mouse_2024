//////////////////////////////////////////////////////
// Code for delivering a sequence of odors to the headfixed mouse
// Cowen and NS&B Students 2023
//////////////////////////////////////////////////////
int n_solenoids = 5;
byte solenoid_pins[] = {9, 10, 11, 12, 13}; // so far, pin/odor port 9 does not work
int port_open_time_ms = 200;

void setup() {
  // initialize digital pin LED_BUILTIN as an output.
  Serial.begin(9600);
  Serial.println("STARTED");
  for (int i = 0;i<n_solenoids;i++){
    pinMode(solenoid_pins[i], OUTPUT);
    digitalWrite(solenoid_pins[i], LOW);
  }
}

// the loop function runs over and over again forever
void loop() {
  for (int i = 0;i<n_solenoids;i++){
    Serial.print("Sol ");
    Serial.println(solenoid_pins[i]);
    digitalWrite(solenoid_pins[i], HIGH);
    delay(port_open_time_ms);                      // wait 
    digitalWrite(solenoid_pins[i], LOW);
    delay(50);                      // wait 
  }
}
