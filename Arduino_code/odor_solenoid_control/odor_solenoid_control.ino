//////////////////////////////////////////////////////
// Code for delivering a sequence of odors to the headfixed mouse
// Cowen and NS&B Students 2023
//////////////////////////////////////////////////////
int n_solenoids = 4;
byte solenoid_pins[] = {10, 11, 12, 13}; 
int sequence_length = 8;
byte odor_sequence[] = {1, 2, 3, 4, 1, 2, 3, 4}; // so far, pin/odor port 9 does not work
byte ttl_odor_code_pin = 5;
int ttl_time_ms = 5;  
//int port_open_time_ms = 1000;
//int port_closed_time_ms = 5000;
int port_open_time_ms = 100;
int port_closed_time_ms = 500;

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
  for (int i = 0;i<sequence_length;i++){
    Serial.print("Sol ");
    Serial.print(i+1);
    Serial.print(", Pin ");
    Serial.println(solenoid_pins[odor_sequence[i] - 1]);
    // Send a code to Neuropixels to indicate the odor type.
    for (int j = 0;j<odor_sequence[i];j++){
      digitalWrite(ttl_odor_code_pin, HIGH);
      delay(ttl_time_ms);
      digitalWrite(ttl_odor_code_pin, LOW);
      delay(ttl_time_ms);
    }
    digitalWrite(solenoid_pins[i], HIGH);
    delay(port_open_time_ms);                      // wait 
    digitalWrite(solenoid_pins[i], LOW);
    delay(port_closed_time_ms);                      // wait 
  }
}
