import cv2
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import math

import os


# Función para calcular el ángulo entre dos vectores
def angle_between(v1, v2):
    dot_product = np.dot(v1, v2)
    norms = np.linalg.norm(v1) * np.linalg.norm(v2)
    cos_theta = dot_product / norms
    return np.arccos(cos_theta)

# Función para calcular la derivada parcial de theta con respecto a x o y
def partial_derivative_theta(v1, v2, var):
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)
    dot_product = np.dot(v1, v2)
    if var == 'x1':
        return (v2[0] * norm_v1 - dot_product * v1[0] / norm_v1) / (norm_v1**2 * norm_v2)
    elif var == 'y1':
        return (v2[1] * norm_v1 - dot_product * v1[1] / norm_v1) / (norm_v1**2 * norm_v2)
    elif var == 'x2':
        return (v1[0] * norm_v2 - dot_product * v2[0] / norm_v2) / (norm_v1 * norm_v2**2)
    elif var == 'y2':
        return (v1[1] * norm_v2 - dot_product * v2[1] / norm_v2) / (norm_v1 * norm_v2**2)


def procesar_img(nombre_imagen):

  vector_center_to_yellow = None
  vector_yellow_to_blue = None

  # Cargar la imagen PNG
  image = cv2.imread(os.path.join(carpeta_imagenes, nombre_imagen))

  # Convertir la imagen a espacio de color HSV
  hsv = cv2.cvtColor(image, cv2.COLOR_BGR2HSV)

  # Definir los límites inferior y superior para el color morado en HSV (ajusta si es necesario)
  lower_purple = np.array([140, 100, 150])
  upper_purple = np.array([160, 255, 255])

  # Crear una máscara que filtre solo los píxeles de color morado
  mask = cv2.inRange(hsv, lower_purple, upper_purple)

  # Aplicar un poco de desenfoque para mejorar la detección de círculos
  mask_blurred = cv2.GaussianBlur(mask, (9, 9), 2)

  # Detectar los círculos usando HoughCircles en la máscara
  circles = cv2.HoughCircles(mask_blurred, cv2.HOUGH_GRADIENT, dp=1.2, minDist=30, param1=50, param2=20, minRadius=7, maxRadius=15)

  # Verificar si se detectaron círculos
  if circles is not None:
      circles = np.round(circles[0, :]).astype("int")

      # Dibujar los círculos detectados en la máscara
      mask_with_circles = cv2.cvtColor(mask, cv2.COLOR_GRAY2BGR)
      for (x, y, r) in circles:
          cv2.circle(mask_with_circles, (x, y), r, (0, 255, 0), 2)  # Dibujar círculos verdes

      # Mostrar la máscara con los círculos detectados
      #plt.imshow(cv2.cvtColor(mask_with_circles, cv2.COLOR_BGR2RGB))
      #plt.title("Círculos detectados en la máscara")
      #plt.axis("off")
      #plt.show()

      # Extraer las coordenadas x, y de los círculos detectados
      x_coords = [x for x, y, r in circles]
      y_coords = [y for x, y, r in circles]

      # Verificar que se hayan detectado al menos 4 círculos
      if len(x_coords) >= 4:

          # Función de error para el ajuste del círculo
          def calc_error(params, x_coords, y_coords):
              x_c, y_c, r = params
              residuals = [(np.sqrt((x - x_c)**2 + (y - y_c)**2) - r)**2 for x, y in zip(x_coords, y_coords)]
              return np.sum(residuals)

        # Estimación inicial mejorada
          x_init = np.mean(x_coords)
          y_init = np.mean(y_coords)
          r_init = np.mean([np.sqrt((x - x_init)**2 + (y - y_init)**2) for x, y in zip(x_coords, y_coords)])
          initial_guess = [x_init, y_init, r_init]

        # Optimización usando mínimos cuadrados
          result = minimize(calc_error, initial_guess, args=(x_coords, y_coords), method='Powell')

          if result.success:
              x_c_opt, y_c_opt, r_opt = result.x
              x_c_opt, y_c_opt, r_opt = int(x_c_opt), int(y_c_opt), int(r_opt)

             # Calcular el error del radio
              errors = [(np.sqrt((x - x_c_opt)**2 + (y - y_c_opt)**2) - r_opt) for x, y in zip(x_coords, y_coords)]
              error_r = np.sqrt(np.sum(np.square(errors))) / len(errors)  # Error promedio del radio

              print(f"Radio ajustado, {r_opt}, px,", end=" ")
              print(f"Error del radio, {error_r:.2f}, px,", end=" ")

              # Dibujar el círculo ajustado en el centro con el radio calculado en color blanco
              cv2.circle(image, (x_c_opt, y_c_opt), r_opt, (255, 255, 255), 2)  # Círculo blanco

              # Detectar el punto amarillo (color ffde59)
              lower_yellow = np.array([25, 100, 200])  # Limite inferior del amarillo en HSV
              upper_yellow = np.array([35, 255, 255])  # Limite superior del amarillo en HSV
              mask_yellow = cv2.inRange(hsv, lower_yellow, upper_yellow)

              # Encontrar contornos del punto amarillo
              contours, _ = cv2.findContours(mask_yellow, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
              if contours:
                  # Obtener el contorno más grande
                  largest_contour = max(contours, key=cv2.contourArea)
                  M = cv2.moments(largest_contour)
                  if M["m00"] > 0:
                      yellow_x = int(M["m10"] / M["m00"])
                      yellow_y = int(M["m01"] / M["m00"])

                      # Dibujar un círculo para marcar el punto amarillo
                      cv2.circle(image, (yellow_x, yellow_y), 15, (255, 222, 89), -1)

                      # Calcular la pendiente de la línea que pasa por el centro del círculo y el punto amarillo
                      dx = yellow_x - x_c_opt
                      dy = yellow_y - y_c_opt

                      # Normalizar la pendiente para extender la línea
                      length = np.sqrt(dx**2 + dy**2)
                      dx /= length
                      dy /= length

                      # Calcular las intersecciones de la línea con el círculo (extendida hacia ambos lados)
                      intersection1_x = int(x_c_opt + r_opt * dx)
                      intersection1_y = int(y_c_opt + r_opt * dy)
                      intersection2_x = int(x_c_opt - r_opt * dx)
                      intersection2_y = int(y_c_opt - r_opt * dy)

                      # Dibujar la línea desde una intersección hasta la otra
                      cv2.line(image, (intersection1_x, intersection1_y), (intersection2_x, intersection2_y), (0, 255, 255), 2)

                      # Detectar el punto azul (color 0319d8)
                      lower_blue = np.array([110, 150, 100])  # Limite inferior del azul en HSV
                      upper_blue = np.array([130, 255, 255])  # Limite superior del azul en HSV
                      mask_blue = cv2.inRange(hsv, lower_blue, upper_blue)

                      # Encontrar contornos del punto azul
                      contours_blue, _ = cv2.findContours(mask_blue, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
                      if contours_blue:
                          # Obtener el contorno más grande
                          largest_contour_blue = max(contours_blue, key=cv2.contourArea)
                          M_blue = cv2.moments(largest_contour_blue)
                          if M_blue["m00"] > 0:
                              blue_x = int(M_blue["m10"] / M_blue["m00"])
                              blue_y = int(M_blue["m01"] / M_blue["m00"])

                              # Dibujar un círculo para marcar el punto azul
                              cv2.circle(image, (blue_x, blue_y), 15, (3, 25, 216), -1)

                              # Dibujar la línea desde el punto amarillo al punto azul
                              cv2.line(image, (yellow_x, yellow_y), (blue_x, blue_y), (255, 0, 0), 2)

                              # Calcular el ángulo entre las dos líneas
                              vector_center_to_yellow = np.array([ x_c_opt - yellow_x,  y_c_opt - yellow_y])
                              vector_yellow_to_blue = np.array([blue_x - yellow_x, blue_y - yellow_y])

                              magnitude_center_to_yellow = np.linalg.norm(vector_center_to_yellow)
                              magnitude_yellow_to_blue = np.linalg.norm(vector_yellow_to_blue)

                              dot_product = np.dot(vector_center_to_yellow, vector_yellow_to_blue)
                              cos_angle = dot_product / (magnitude_center_to_yellow * magnitude_yellow_to_blue)
                              cos_angle = np.clip(cos_angle, -1.0, 1.0)

                              angle_radians = np.arccos(cos_angle)
                              angle_degrees = np.degrees(angle_radians)

                              print(f"El ángulo entre las dos líneas es, {angle_degrees:.2f}, grados,", end=" ")

                              # Dibujar vectores
                              arrow_end_center_to_yellow = (x_c_opt, y_c_opt)
                              cv2.arrowedLine(image, (yellow_x, yellow_y), arrow_end_center_to_yellow, (0, 255, 255), 2, tipLength=0.05)

                              arrow_end_yellow_to_blue = (blue_x, blue_y)
                              cv2.arrowedLine(image, (yellow_x, yellow_y), arrow_end_yellow_to_blue, (255, 0, 0), 2, tipLength=0.05)

                          else:
                              print("No se pudo encontrar el punto azul.")
                      else:
                          print("No se detectó el punto azul.")
                  else:
                      print("No se pudo encontrar el punto amarillo.")
              else:
                  print("No se detectó el punto amarillo.")
          else:
              print("No se pudo ajustar el círculo. Razón:", result.message)
      else:
          print("No se detectaron al menos 4 círculos.")
  else:
      print("No se detectaron círculos.")

  # Error de 0.5 px para cada componente
  error_px = 0.5

  # Cálculo del ángulo entre los dos vectores
  theta = angle_between(vector_center_to_yellow, vector_yellow_to_blue)

  # Calcular las derivadas parciales de theta con respecto a las variables
  dtheta_dx1 = partial_derivative_theta(vector_center_to_yellow, vector_yellow_to_blue, 'x1')
  dtheta_dy1 = partial_derivative_theta(vector_center_to_yellow, vector_yellow_to_blue, 'y1')
  dtheta_dx2 = partial_derivative_theta(vector_center_to_yellow, vector_yellow_to_blue, 'x2')
  dtheta_dy2 = partial_derivative_theta(vector_center_to_yellow, vector_yellow_to_blue, 'y2')

  # Aplicar la fórmula de propagación de errores para obtener la incertidumbre del ángulo
  uncertainty_theta = np.sqrt(
      (dtheta_dx1 * error_px)**2 +
      (dtheta_dy1 * error_px)**2 +
      (dtheta_dx2 * error_px)**2 +
      (dtheta_dy2 * error_px)**2
  )

  print(f"Ángulo calculado, {theta}, rad,", end=" ")
  print(f"Incertidumbre del ángulo, {uncertainty_theta}, rad")

  # Mostrar la imagen original con los círculos y vectores dibujados
  #plt.figure(figsize=(12, 8), dpi=300)
  #plt.imshow(cv2.cvtColor(image, cv2.COLOR_BGR2RGB))
  #plt.title(f"Imagen con los círculos y vectores trazados {nombre_imagen}")
  #plt.axis("off")
  #plt.savefig(f'/content/res/vec_{nombre_imagen}.png', dpi=300)
  #print(f"Guardada la imagen: vec_{nombre_imagen}.png")
  #plt.show()

# Procesar cada imagen en la carpeta
for i in range(1, 7):  # Cambiar a range(1, 7) para iterar de 1 a 6

    # Ruta de la carpeta donde se encuentran las imágenes
    carpeta_imagenes = f"/content/{i}"  # Cambia esto a la ruta correcta

    # Obtener una lista de todas las imágenes en la carpeta
    imagenes = [f for f in os.listdir(carpeta_imagenes) if f.endswith('.png')]

    print(f"Medición: {i}")

    for nombre_imagen in imagenes:
        print(f"foto, {nombre_imagen},", end=" ")
        procesar_img(nombre_imagen)
    print("")