# Spinning Top Java Project

This project contains a simple spinning top simulation using JOGL.

## Building

The project uses Maven for building. To compile and package the application run:

```bash
mvn package
```

This will produce a jar file in the `target` directory with the required manifest to run the viewer.

## Running

Run the jar with:

```bash
java -jar target/spinning-top-1.0-SNAPSHOT.jar
```

Make sure JOGL native libraries are available on your system when running the program.
