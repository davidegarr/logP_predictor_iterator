import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import tensorflow as tf
import matplotlib.pyplot as plt 
from sklearn.metrics import r2_score
from joblib import dump

def main():
    #load the data
    data = pd.read_csv("clean_data.csv")

    #features and labels
    features = ["MolWt","NumHAcceptors", "NumHDonors","NumHeteroatoms","NumRotatableBonds","NumAromaticCarbocycles","NumAromaticHeterocycles","NumSaturatedHeterocycles","NumSaturatedCarbocycles","FractionCSP3","TPSA"]
    X = data[features]
    y = data["AlogP"]

    #split the data into training and validation 
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    model = tf.keras.Sequential([
    tf.keras.layers.Dense(256, activation='relu', input_shape=(X_train_scaled.shape[1],)),
    tf.keras.layers.Dense(128, activation='relu'),
    tf.keras.layers.Dense(64, activation='relu'),
    tf.keras.layers.Dense(32, activation='relu'),
    tf.keras.layers.Dense(1)
    ])

    model.compile(optimizer='adam', loss='mse', metrics=['mae'])

    # Fit the model
    history = model.fit(X_train_scaled, y_train, epochs=25, batch_size=512, validation_split=0.2)

    # Evaluate the model
    test_loss, test_mae = model.evaluate(X_test_scaled, y_test)
    print(f"Test Loss: {test_loss}, Test MAE: {test_mae}")

    # Generate predictions on test set
    y_pred = model.predict(X_test_scaled).flatten()

    #SAve the model
    model.save("./predictor_v1")
    dump(scaler, "scaler.joblib")
 

    # Compute R^2 value
    r2 = r2_score(y_test, y_pred)

    # Plot actual vs predicted values
    plt.scatter(y_test, y_pred)
    plt.xlabel("Actual AlogP Values")
    plt.ylabel("Predicted AlogP Values")
    plt.title(f"Actual vs Predicted AlogP Values (R^2: {r2:.2f})")

    # Draw a diagonal line for perfect predictions
    min_val = min(min(y_test), min(y_pred))
    max_val = max(max(y_test), max(y_pred))
    plt.plot([min_val, max_val], [min_val, max_val], 'r--')

    # Display R^2 value on the plot
    plt.text(min_val, max_val, f'R^2 = {r2:.2f}', horizontalalignment='left', verticalalignment='bottom')

    plt.show()


main()