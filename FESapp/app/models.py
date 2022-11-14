from django.db import models

# Create your models here.
class Item(models.Model):
    """ a model with name, descriptiom, price, image and time along with status """
    name = models.CharField(max_length=100)
    description = models.TextField()
    price = models.IntegerField()
    image = models.ImageField(upload_to='pics')
    time = models.DateTimeField(auto_now_add=True)
    status = models.BooleanField(default=True)

    def __str__(self):
        return self.name

    #check if a number is prime
    def is_prime(self):
        if self.price > 1:
            for i in range(2, self.price):
                if (self.price % i) == 0:
                    return False
            else:
                return True
        else:
            return False


    def add(self, num1, num2):
        return num1 + num2

    def sub(self, num1, num2):
        return num1 - num2
    
    def mul(self, num1, num2):  
        return num1 * num2

    def div(self, num1, num2):
        return num1 / num2
    
    def mod(self, num1, num2):
        return num1 % num2

    mergesort = lambda self, arr: arr if len(arr) <= 1 else self.merge(self.mergesort(arr[:len(arr)//2]), self.mergesort(arr[len(arr)//2:]))
    
    
    