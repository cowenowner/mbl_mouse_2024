int n_solenoids = 5;
byte solenoid_pins[] = {9, 10, 11, 12, 13};

void setup() {
  // initialize digital pin LED_BUILTIN as an output.
  Serial.print("STARTED");
  for (int i = 0;i<n_solenoids;i++){
    pinMode(solenoid_pins[i], OUTPUT);
    digitalWrite(solenoid_pins[i], LOW);

  }
}

// the loop function runs over and over again forever
void loop() {

  for (int i = 0;i<n_solenoids;i++){
    Serial.print("Sol ");
    Serial.print(solenoid_pins[i]);
    Serial.print(",");
    digitalWrite(solenoid_pins[i], HIGH);
    delay(1000);                      // wait for a second
    digitalWrite(solenoid_pins[i], LOW);
    delay(50);                      // wait for a second

  }
}
