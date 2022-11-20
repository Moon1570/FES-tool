from django.shortcuts import render
from django.http import HttpResponse
from django.http import JsonResponse


# Create your views here.
def home(request):
    return render(request, 'home.html', {'name':'Moon'})

def calc(request):
    #get name from request
    name = request.POST.get('name')
    weight = int(request.POST.get('weight'))
    print(weight)
    intervalTime = int(request.POST.get('intervalTime'))

    #print to console
    print(name)
    print(weight)
    print(intervalTime)



    return JsonResponse({'name':name, 'weight':weight, 'intervalTime':intervalTime})
